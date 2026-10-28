#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PROTENIX_PREDICT {
    tag { "${dataset_name}_${chunk_file.simpleName}" }

    label 'gpu'

    input:
    tuple val(dataset_name), path(chunk_file)

    output:
    tuple val(dataset_name), path("protenix_out/${dataset_name}"), emit: protenix_out

    script:
    """
    # module load cuda/13.0.2
    # module load python/3.10
    # source ~/miniconda3/etc/profile.d/conda.sh
    # conda activate protenix_env
    export PROTENIX_CACHE=${params.protenix_cache}
    export PROTENIX_ROOT_DIR=${params.protenix_cache}
    export MPLCONFIGDIR=\$(mktemp -d)
    export TORCH_EXTENSIONS_DIR=\$PWD/torch_extensions
    export XDG_CACHE_HOME=\$PWD/.cache
    export TMPDIR=\$(mktemp -d)
    export MPLCONFIGDIR=\$PWD/.config

    mkdir -p "\$TORCH_EXTENSIONS_DIR" "\$XDG_CACHE_HOME"  "\$MPLCONFIGDIR"

    export HOME=\$PWD/home

    # This is too hard coded and should be updated to be more flexible.

    mkdir -p protenix_out/${dataset_name}
    
    echo \$PROTENIX_ROOT_DIR
   
    python /app/runner/batch_inference.py \
        --input ${chunk_file} \
        --out_dir protenix_out/${dataset_name} \
        --seeds 101 \
        --model_name "protenix_mini_esm_v0.5.0" \
        --use_msa false \
        --sample 1
    """
}
