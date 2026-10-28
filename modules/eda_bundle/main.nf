process EDA_BUNDLE {
    tag "eda_bundle"
    // Fixed name: overwritten each run, like everything else in output_dir.
    // eda/run_config.* and eda/biobank_config.* are also published loose,
    // next to the report.
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    // Rebuild the published eda/ layout so the report's profiles/ links resolve
    path report_files, stageAs: 'eda/*'
    path profiles,     stageAs: 'eda/profiles/*'
    val  run_info
    path biobank_config, stageAs: 'input/biobank_config.yaml'
    path eda_config,     stageAs: 'input/eda_config.json'

    output:
    path "eda_report.zip",       emit: zip
    path "eda/run_config.html",  emit: run_config_html
    path "eda/run_config.json",  emit: run_config_json
    path "eda/biobank_config.*", emit: biobank_config

    script:
    """
    cat > input/run_info.json <<'RUN_INFO_EOF'
${run_info}
RUN_INFO_EOF
    run_config.py \\
        --run_info       input/run_info.json \\
        --biobank_config input/biobank_config.yaml \\
        --eda_config     input/eda_config.json \\
        --output_dir     eda
    python3 -m zipfile -c eda_report.zip eda/
    """

    stub:
    """
    touch eda_report.zip eda/run_config.html eda/run_config.json \\
          eda/biobank_config.yaml eda/biobank_config.html
    """
}
