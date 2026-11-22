// OpenDNA as Nextflow processes.
//
// The point of running OpenDNA under a workflow engine is parallelism across
// a few hundred candidate sequences, and the thing that makes it safe is that
// the CLI already answers the questions a scheduler cannot: which models are
// installed, whether they verify, what each failure means, and whether two
// parts were produced the same way.
//
// So these processes are thin. They do not re-implement any of that; they
// wire it up.
//
// Three things are deliberate.
//
// `OPENDNA_CHECK` runs once, before anything is scheduled, and fails the whole
// run if the installed heads are not the ones in the lockfile. Discovering
// that after four hundred tasks have run is discovering it too late.
//
// `errorStrategy` distinguishes the exit codes. 3 means the installation is
// wrong, 4 means the data is wrong, 5 means the install does not reproduce its
// own recorded results: retrying any of those just burns the queue. Only 1,
// which is an unanticipated failure and might be a full disk or an evicted
// node, is retried.
//
// `OPENDNA_MERGE` uses `opendna merge` rather than `cat`, because it reads
// each part's run manifest and refuses to join parts produced on different
// devices or with a different head build. On a cluster where nodes are
// upgraded one at a time, that is a real thing to catch.

nextflow.enable.dsl = 2

// Verify the installation before any work is scheduled.
process OPENDNA_CHECK {
    label 'opendna_small'

    input:
    path models
    path lockfile

    output:
    path 'models.ok', emit: ok

    script:
    """
    opendna --models ${models} models --check ${lockfile}
    opendna --models ${models} selftest
    touch models.ok
    """
}

// Cut the input into parts balanced by residue count.
process OPENDNA_SPLIT {
    label 'opendna_small'

    input:
    path fasta
    val parts

    output:
    path 'parts/part-*.fasta', emit: parts

    script:
    """
    opendna split -i ${fasta} --out-dir parts --parts ${parts}
    """
}

// Predict one part.
//
// `--offline` is on: a staged models directory holds the embedding model, and
// a task that reaches for the network on a compute node either stalls or
// downloads 568 MB once per task. If the bundle is incomplete this fails
// immediately with exit 3 instead.
process OPENDNA_PREDICT {
    label 'opendna_predict'
    // Only an unanticipated failure is worth retrying. 3, 4 and 5 are
    // conditions that will hold just as firmly the second time.
    errorStrategy { task.exitStatus == 1 && task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path part
    path models
    val task_name
    val batch_tokens

    output:
    path "${part.baseName}.csv", emit: results
    path "${part.baseName}.manifest.json", emit: manifest

    script:
    """
    opendna --models ${models} --offline --batch-tokens ${batch_tokens} \\
        predict --task ${task_name} -i ${part} -o ${part.baseName}.csv
    """
}

// Join the parts, refusing a mixture.
process OPENDNA_MERGE {
    label 'opendna_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path results
    path manifests
    val name

    output:
    path "${name}.csv", emit: merged

    script:
    """
    opendna merge --output ${name}.csv \$(ls part-*.csv | sort)
    """
}

workflow OPENDNA {
    take:
    fasta       // path to the input FASTA
    models      // path to a staged models directory
    lockfile    // path to models.lock
    task_name   // 'disorder' | 'burial' | 'solubility'
    parts       // how many parts to split into

    main:
    OPENDNA_CHECK(models, lockfile)

    // Nothing is scheduled until the check has passed.
    OPENDNA_SPLIT(fasta, parts)

    OPENDNA_PREDICT(
        OPENDNA_SPLIT.out.parts.flatten(),
        models,
        task_name,
        params.batch_tokens
    )

    OPENDNA_MERGE(
        OPENDNA_PREDICT.out.results.collect(),
        OPENDNA_PREDICT.out.manifest.collect(),
        task_name
    )

    emit:
    merged = OPENDNA_MERGE.out.merged
}
