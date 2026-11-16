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
include { EDA_STATS }          from './modules/eda_stats/main.nf'
include { EDA_REPORT }         from './modules/eda_report/main.nf'
include { EDA_PROFILE }        from './modules/eda_profile/main.nf'
include { EDA_BUNDLE }         from './modules/eda_bundle/main.nf'
include { FILTER_PHENOTYPES }  from './modules/filter_phenotypes/main.nf'

// Advanced phenotypes with config warnings (manifest `warnings` column);
// counted from ADVANCED_CODED.out.warnings, reported in workflow.onComplete
advancedCodedWarned = []

workflow {
    // -------------------------------------------------------------------------
    // Run summary: logged here and, with --run_eda, written to
    // eda/run_config.html by EDA_BUNDLE. Keep the two in sync by editing only
    // this map.
    // -------------------------------------------------------------------------
    def git_rev = workflow.commitId ?: {
        try {
            def rev = ['git', '-C', "${projectDir}", 'rev-parse', '--short', 'HEAD'].execute()
            rev.waitFor()
            if (rev.exitValue() != 0) return null
            def dirty = ['git', '-C', "${projectDir}", 'status', '--porcelain'].execute().text.trim()
            return rev.text.trim() + (dirty ? ' (uncommitted changes)' : '')
        } catch (Exception e) {
            return null
        }
    }()
    def run_info = [
        'RUN': [
            'run as'           : workflow.commandLine,
            'run location'     : "${launchDir}",
            'started at'       : "${workflow.start}",
            'profile'          : workflow.profile,
            'resumed'          : workflow.resume,
            'session id'       : "${workflow.sessionId}",
            'nextflow version' : "${nextflow.version}",
            'pipeline revision': git_rev ?: 'unknown',
            'container'        : workflow.container ? "${workflow.container}" : 'none',
            'work dir'         : "${workflow.workDir}",
        ],
        'WORKFLOW PARAMETERS': [
            '--biobank_config'         : params.biobank_config ?: 'NOT SET',
            '--run_gather'             : params.run_gather ? 'ENABLED' : 'DISABLED',
            '--gather_format'          : params.gather_format,
            '--run_risk_classification': params.run_risk_classification ? 'ENABLED' : 'DISABLED',
            '--run_eda'                : params.run_eda ? 'ENABLED' : 'DISABLED',
            '--filter_config'          : params.filter_config ?: 'NOT SET (no filtering)',
            '--event_mask_config'      : params.event_mask_config ?: 'NOT SET (no masking)',
            '--sample_list'            : params.sample_list ?: 'NOT SET (all samples)',
            '--coded_chunk_size'       : params.coded_chunk_size ?: 0,
            '--chunk_size'             : params.chunk_size ?: 0,
            '--eda_batch_size'         : params.eda_batch_size,
        ],
        'OUTPUT DIRECTORIES': [
            '--output_dir': params.output_dir ?: 'NOT SET',
        ],
    ]
    def log_lines = ['  NEXTFLOW - DSL2 - PHENOTYPING - P I P E L I N E', '  ' + '=' * 50]
    run_info.each { section, rows ->
        if (section != 'RUN') log_lines += ['', "  ${section}:", '  ' + '=' * 50]
        rows.each { k, v -> log_lines << String.format('  %-25s : %s', k, v) }
    }
    log.info((log_lines + ['  ' + '=' * 50]).join('\n'))

    // Removed from nextflow.config: they never had any effect (phenotype
    // settings live in the biobank YAML). Warn instead of silently ignoring.
    def removed_params = [
        'run_preprocess', 'run_quantitative', 'run_categorical', 'input', 'long_format_input',
        'quantitative_params', 'categorical_params', 'input_format', 'data_type',
        'sample_id_col', 'phenotype_col', 'value_col', 'date_col', 'reference_date',
        'outlier_method', 'outlier_iqr_multiplier', 'outlier_zscore_sd', 'outlier_mode',
        'outlier_phenotypes', 'scale_phenotypes', 'scale_method', 'quant_stats',
        'quant_phenotypes', 'cat_phenotypes', 'missing_as_control', 'min_occurrences',
    ].findAll { params.containsKey(it) }
    if (removed_params) {
        log.warn "Ignoring unused parameter(s): ${removed_params.collect { '--' + it }.join(', ')}. " +
                 "Phenotype settings belong in the biobank YAML (--biobank_config)."
    }

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
    if (!params.sample_list) {
        log.warn "--sample_list NOT SET: coded/advanced_coded phenotypes have no full sample roster, " +
                 "so missing_as_control cannot write 0s (non-cases will be NA)."
    }
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
    BASIC_CODED(diag_pairs, ch_sample_list)

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

    ADVANCED_CODED(advanced_inputs, ch_sample_list)

    // Unmatched-pattern / 0-case warnings go in the manifest's `warnings`
    // column; here we only count the phenotypes that have any, for a one-line
    // pointer at the end of the run. Read from task outputs, so cached tasks
    // count again on -resume.
    ADVANCED_CODED.out.warnings
        .filter { f -> f.text.trim() }
        .subscribe { f -> advancedCodedWarned << f.name.replace('.warnings.txt', '') }

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
    // Phenotype manifest — one row per output column of every phenotype
    //   (column_name, phenotype, data_type, code) → ${output_dir}/phenotype_manifest.tsv
    //   Sorted by data_type, then phenotype, then column; rows within one
    //   phenotype keep their output order (e.g. BMI_mean before BMI_std).
    //   Risk classification columns are not included (separate risk_matrix.tsv).
    //   With --run_eda, EDA_REPORT publishes this manifest instead, with an
    //   added `flags` column, so it is not stored here (no double write).
    // -------------------------------------------------------------------------
    def manifest_opts = [name: 'phenotype_manifest.tsv', keepHeader: true, skip: 1,
                         sort: { f ->
                             def lines = f.readLines()
                             if (lines.size() < 2) return ''
                             def row = lines[1].split('\t', -1)
                             [row[2], row[1], row[0]].join('\t')
                         }]
    if (!params.run_eda) manifest_opts.storeDir = params.output_dir
    manifest_ch = BASIC_QUANTITATIVE.out.manifest
        .mix(BASIC_CATEGORICAL.out.manifest)
        .mix(BASIC_CODED.out.manifest)
        .mix(ADVANCED_CODED.out.manifest)
        .collectFile(manifest_opts)

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

    // -------------------------------------------------------------------------
    // Step 6: Optional EDA / QC report
    //   EDA_STATS  : per-column metrics for every result, in batches, against
    //                the covariates named in the YAML's eda: block.
    //   EDA_REPORT : flags, triage HTML, restriction suggestions, and the
    //                bounded profile selection (eda cascade + top-k flagged).
    //   EDA_PROFILE: one ydata report per selected group, in parallel.
    //   EDA_BUNDLE : writes eda/run_config.html (run parameters + resolved
    //                YAML settings), then zips eda/ to {output_dir}/eda_report.zip.
    //   Reads per-phenotype results directly — does not need --run_gather.
    // -------------------------------------------------------------------------
    if (params.run_eda) {
        def eda_yaml = new org.yaml.snakeyaml.Yaml()
            .load(new File(params.biobank_config).text)
        def covariate_names = ((eda_yaml?.eda?.covariates ?: []) as List)
            .collect { it.toString() }.toSet()

        eda_results = quant_results.mix(cat_results).mix(diag_results)

        // Result files are named <phenotype>.quant.tsv / <phenotype>.cat.tsv
        covariate_files = quant_results.mix(cat_results)
            .filter { f -> covariate_names.contains(f.name.replaceAll(/\.(quant|cat)\.tsv$/, '')) }
            .toSortedList { a, b -> a.name <=> b.name }

        eda_config_ch = RESOLVE_CONFIG.out.eda_config.first()

        // Sort before batching: arrival order varies run to run, and unstable
        // batches would miss the cache on -resume (cascading to EDA_REPORT)
        EDA_STATS(
            eda_results
                .toSortedList { a, b -> a.name <=> b.name }
                .flatMap()
                .collate(params.eda_batch_size ?: 200),
            covariate_files,
            eda_config_ch,
            ch_sample_list
        )

        eda_stats_ch = EDA_STATS.out.stats
            .collectFile(name: 'eda_stats.tsv', keepHeader: true, skip: 1)

        EDA_REPORT(
            eda_stats_ch,
            manifest_ch,
            eda_config_ch,
            // All named eda_cases.batch.tsv; sort by work-dir path for a stable cache key
            EDA_STATS.out.cases.toSortedList { a, b -> a.toString() <=> b.toString() },
            // Identical in every batch; take a fixed one so the cache key is stable
            EDA_STATS.out.covariates.toSortedList { a, b -> a.toString() <=> b.toString() }.map { it[0] }
        )

        // Value channel: read by both the splitCsv below and every EDA_PROFILE task
        profile_selection_ch = EDA_REPORT.out.profile_selection.first()

        // Join the selection back to result files by name, one task per group
        results_by_name = eda_results.map { f -> [f.name, f] }
        profile_inputs = profile_selection_ch
            .splitCsv(header: true, sep: '\t')
            .map { r -> [r.source_file, r.group] }
            .unique()
            .combine(results_by_name, by: 0)
            .map { name, group, f -> [group, f] }
            .groupTuple()

        EDA_PROFILE(
            profile_inputs,
            profile_selection_ch,
            eda_config_ch
        )

        // One zip of the eda/ folder for end users. collect waits for the report;
        // toSortedList emits [] when EDA_PROFILE produced nothing (profiler
        // missing), so the bundle still runs without profiles
        EDA_BUNDLE(
            EDA_REPORT.out.summary
                .mix(EDA_REPORT.out.flags, EDA_REPORT.out.html,
                     EDA_REPORT.out.profile_selection, EDA_REPORT.out.restrictions,
                     EDA_REPORT.out.filter_suggested, EDA_REPORT.out.covariates,
                     EDA_REPORT.out.pairwise)
                .collect(sort: { a, b -> a.name <=> b.name }),
            EDA_PROFILE.out.html
                .toSortedList { a, b -> a.name <=> b.name },
            // Changes every run (start time), so only this cheap step reruns on -resume
            Channel.value(groovy.json.JsonOutput.toJson(
                run_info + [start_epoch_ms: workflow.start.toInstant().toEpochMilli()])),
            file(params.biobank_config),
            eda_config_ch
        )
    }

    // -------------------------------------------------------------------------
    // Step 7: Optional filter (--filter_config)
    //   Post-processing only: reads phenotype_table.wide.tsv and the EDA
    //   tables, writes filtered/phenotype_table.filtered.wide.tsv + log.
    //   Raw outputs are never modified. Start from eda/filter_suggested.yaml.
    // -------------------------------------------------------------------------
    if (params.filter_config) {
        if (!params.run_gather || !params.run_eda) {
            error "ERROR: --filter_config requires --run_gather and --run_eda."
        }
        if ((params.gather_format ?: 'wide') == 'long') {
            error "ERROR: --filter_config needs the wide table (--gather_format wide|both)."
        }
        FILTER_PHENOTYPES(
            GATHER.out.wide,
            Channel.value(file(params.filter_config, checkIfExists: true)),
            EDA_REPORT.out.covariates,
            EDA_REPORT.out.summary,
            EDA_REPORT.out.manifest
        )
    }
}

// One-line pointer at the end of the run, where it won't scroll past
workflow.onComplete {
    if (advancedCodedWarned) {
        log.warn "${advancedCodedWarned.size()} advanced_coded phenotype(s) have config warnings " +
                 "— see the `warnings` column in ${params.output_dir}/phenotype_manifest.tsv"
    }
}
