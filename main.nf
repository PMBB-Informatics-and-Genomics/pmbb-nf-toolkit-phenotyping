#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { RESOLVE_CONFIG }    from './modules/resolve_config/main.nf'
include { PREPROCESS as PREPROCESS_DIAG }  from './modules/preprocess/main.nf'
include { PREPROCESS as PREPROCESS_OTHER } from './modules/preprocess/main.nf'
include { BASIC_CODED }    from './modules/basic_coded/main.nf'
include { ADVANCED_CODED } from './modules/advanced_coded/main.nf'
include { BASIC_QUANTITATIVE } from './modules/basic_quantitative/main.nf'
include { BASIC_CATEGORICAL }  from './modules/basic_categorical/main.nf'
include { GATHER }             from './modules/gather/main.nf'
include { EVENT_MASK }         from './modules/event_mask/main.nf'
include { APPLY_MASK }         from './modules/apply_mask/main.nf'
include { RISK_CLASSIFICATION } from './modules/risk_classification/main.nf'

workflow {
    // -------------------------------------------------------------------------
    // Log the parameters for debugging
    log.info([
        "  NEXTFLOW - DSL2 - PHENOTYPING - P I P E L I N E",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "run as", workflow.commandLine),
        String.format("  %-25s : %s", "run location", launchDir),
        String.format("  %-25s : %s", "started at", workflow.start),
        "",
        "  WORKFLOW PARAMETERS:",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "--biobank_config", params.biobank_config ?: "NOT SET"),
        String.format("  %-25s : %s", "--run_gather", params.run_gather ? "ENABLED" : "DISABLED"),
        String.format("  %-25s : %s", "--run_risk_classification", params.run_risk_classification ? "ENABLED" : "DISABLED"),
        String.format("  %-25s : %s", "--event_mask_config", params.event_mask_config ?: "NOT SET (no masking)"),
        String.format("  %-25s : %s", "--run_quantitative", params.run_quantitative ? "ENABLED" : "DISABLED"),
        String.format("  %-25s : %s", "--run_categorical", params.run_categorical ? "ENABLED" : "DISABLED"),
        String.format("  %-25s : %s", "--coded_chunk_size", params.coded_chunk_size ?: 0),
        String.format("  %-25s : %s", "--chunk_size", params.chunk_size ?: 0),
        String.format("  %-25s : %s", "--sample_list", params.sample_list ?: "NOT SET (all samples)"),
        "",
        "  OUTPUT DIRECTORIES:",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "--output_dir", params.output_dir ?: "NOT SET"),
            "  " + "=" * 50
    ].join("\n"))

    // -------------------------------------------------------------------------
    // Step 1: Resolve biobank config → per-phenotype JSON files
    // -------------------------------------------------------------------------
    if (!params.biobank_config) {
        error "ERROR: --biobank_config is required."
    }
    config_ch = Channel.fromPath(params.biobank_config, checkIfExists: true)
    ch_sample_list = params.sample_list
        ? Channel.value(file(params.sample_list, checkIfExists: true))
        : Channel.value([])
    RESOLVE_CONFIG(config_ch)

    // -------------------------------------------------------------------------
    // Step 2: Preprocess each phenotype in parallel
    //   RESOLVE_CONFIG emits separate channels per data_type type.
    //   PREPROCESS takes one JSON and emits: tuple(json, [*.long.tsv])
    //   .transpose() expands into individual (json, tsv) pairs.
    // -------------------------------------------------------------------------
    // Diagnostic configs get their own process alias for higher resource allocation
    coded_jsons  = RESOLVE_CONFIG.out.coded_configs.flatten()
    other_jsons = RESOLVE_CONFIG.out.quant_configs
        .mix(RESOLVE_CONFIG.out.cat_configs)
        .flatten()

    PREPROCESS_DIAG(coded_jsons, ch_sample_list)
    PREPROCESS_OTHER(other_jsons, ch_sample_list)

    // Expand multi-file coded outputs into individual (json, tsv) pairs
    preprocessed_pairs = PREPROCESS_DIAG.out.with_config
        .mix(PREPROCESS_OTHER.out.with_config)
        .transpose()

    // Tag each basic coded long.tsv with its phenotype_name
    // (defined here so it is available to the event_mask block below)
    basic_diag_tagged = PREPROCESS_DIAG.out.with_config
        .flatMap { cfg_json, tsvs ->
            def pheno_name = new groovy.json.JsonSlurper().parseText(cfg_json.text).phenotype_name
            (tsvs instanceof List ? tsvs : [tsvs]).collect { tsv -> [pheno_name, tsv] }
        }

    // -------------------------------------------------------------------------
    // Step 2b: Optional event masking
    //   EVENT_MASK runs once per event source phenotype, reads raw long.tsvs,
    //   computes mask windows from event ICD codes. APPLY_MASK filters
    //   long.tsvs for phenotypes that opt in via event_mask: [...] in config.
    //   Raw-data window computation: EVENT_MASK always reads unmasked data
    //   (Solution 1). Self-masking (source phenotype also masked) is valid
    //   and handled gracefully — see design docs.
    //
    //   Phenotypes WITHOUT a real date_col are NOT masked (apply_mask skips
    //   them and logs a warning). All masking uses real observation dates.
    // -------------------------------------------------------------------------
    if (params.event_mask_config) {
        event_mask_config_ch = Channel.fromPath(params.event_mask_config, checkIfExists: true)

        // Parse source phenotype names from event_mask config
        def event_mask_yaml = new org.yaml.snakeyaml.Yaml()
            .load(new File(params.event_mask_config).text)
        def event_source_names = ((event_mask_yaml?.events?.values()
            ?.collect { it.source } ?: []) as List).toSet()

        // Collect source long.tsvs grouped by source phenotype name
        event_source_inputs = basic_diag_tagged
            .filter { pheno_name, tsv -> event_source_names.contains(pheno_name) }
            .groupTuple()   // [source_name, [tsv1, tsv2, ...]]

        EVENT_MASK(
            event_source_inputs,
            event_mask_config_ch
        )

        // Merge per-source interval files into one mask_intervals.tsv
        mask_intervals_ch = EVENT_MASK.out.intervals
            .collectFile(name: 'mask_intervals.tsv', keepHeader: true, skip: 1)

        // Split preprocessed pairs: those needing masking vs those that don't
        preprocessed_pairs
            .map { cfg_json, tsv ->
                def cfg = new groovy.json.JsonSlurper().parseText(cfg_json.text)
                def needs_mask = (cfg.event_mask ?: []).size() > 0
                [needs_mask, cfg.data_type, cfg_json, tsv]
            }
            .branch {
                masked:   it[0] == true
                unmasked: true
            }
            .set { mask_routing }

        APPLY_MASK(
            mask_routing.masked.map { nm, data_type, json, tsv -> [json, tsv] },
            mask_intervals_ch
        )

        masked_data_type_pairs = APPLY_MASK.out.with_config
            .map { cfg_json, tsv ->
                def data_type = new groovy.json.JsonSlurper().parseText(cfg_json.text).data_type
                [data_type, cfg_json, tsv]
            }

        unmasked_data_type_pairs = mask_routing.unmasked
            .map { nm, data_type, json, tsv -> [data_type, json, tsv] }

        all_data_type_pairs = masked_data_type_pairs.mix(unmasked_data_type_pairs)

    } else {
        all_data_type_pairs = preprocessed_pairs
            .map { cfg_json, tsv ->
                def data_type = new groovy.json.JsonSlurper().parseText(cfg_json.text).data_type
                [data_type, cfg_json, tsv]
            }
    }

    all_data_type_pairs
        .branch {
            quant: it[0] == 'quantitative'
            cat:   it[0] == 'categorical'
            diag:  it[0] == 'coded'
        }
        .set { routed }

    quant_pairs = routed.quant.map { c, j, t -> [j, t] }
    cat_pairs   = routed.cat.map   { c, j, t -> [j, t] }
    diag_pairs  = routed.diag.map  { c, j, t -> [j, t] }

    // -------------------------------------------------------------------------
    // Step 3a: Diagnostic codes
    // -------------------------------------------------------------------------
    BASIC_CODED(diag_pairs)

    // -------------------------------------------------------------------------
    // Step 3b: Advanced coded phenotypes
    //   Each advanced phenotype references one or more basic coded
    //   phenotypes by name (from_phenotype). We tag each basic long.tsv with
    //   its parent phenotype_name, then cross with advanced configs to collect
    //   the relevant files for each advanced phenotype.
    // -------------------------------------------------------------------------
    // Cross each advanced config with all basic tagged outputs, keep matching ones,
    // then group all matching tsvs per advanced config
    advanced_inputs = RESOLVE_CONFIG.out.advanced_coded_configs
        .flatten()
        .combine(basic_diag_tagged)
        .filter { adv_json, pheno_name, tsv ->
            def sources = new groovy.json.JsonSlurper().parseText(adv_json.text).sources
            sources != null && sources.any { k, v -> v.from_phenotype == pheno_name }
        }
        .map { adv_json, pheno_name, tsv -> [adv_json.toString(), tsv] }
        .groupTuple()
        .map { adv_json_str, tsvs -> [file(adv_json_str), tsvs.flatten()] }

    ADVANCED_CODED(advanced_inputs)

    diag_results = BASIC_CODED.out.result
        .mix(ADVANCED_CODED.out.result)

    // -------------------------------------------------------------------------
    // Step 3c: Quantitative and categorical basic phenotyping
    // -------------------------------------------------------------------------
    BASIC_QUANTITATIVE(quant_pairs)
    quant_results = BASIC_QUANTITATIVE.out.result

    BASIC_CATEGORICAL(cat_pairs)
    cat_results = BASIC_CATEGORICAL.out.result

    // -------------------------------------------------------------------------
    // Step 4: Optional gather — outer-join all per-phenotype TSVs
    // -------------------------------------------------------------------------
    if (params.run_gather) {
        all_results = quant_results
            .mix(cat_results)
            .mix(diag_results)
            .collect()
        GATHER(all_results)
    }

    // -------------------------------------------------------------------------
    // Step 5: Optional risk classification
    //   Reads gathered wide TSV + resolved risk config from RESOLVE_CONFIG.
    //   Config is defined in pmbb_config.yaml under data_types.quantitative
    //   and/or per phenotype.
    // -------------------------------------------------------------------------
    if (params.run_risk_classification) {
        if (!params.run_gather) {
            error "ERROR: --run_risk_classification requires --run_gather (needs wide output)."
        }
        RISK_CLASSIFICATION(GATHER.out.wide, RESOLVE_CONFIG.out.risk_config)
    }
}
