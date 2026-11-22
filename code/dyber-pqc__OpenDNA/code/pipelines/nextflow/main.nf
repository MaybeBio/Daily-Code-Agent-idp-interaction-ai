// A runnable example of the OpenDNA module.
//
//   nextflow run main.nf \
//       --fasta candidates.fasta \
//       --models /mnt/shared/opendna-models \
//       --lockfile /mnt/shared/opendna-models/models.lock \
//       --task disorder \
//       --parts 20
//
// The models directory is one produced by `opendna stage`, which holds the
// heads, the embedding model, the lockfile and the checksums. Put it on shared
// storage the compute nodes can read; nothing in the run touches the network.

nextflow.enable.dsl = 2

include { OPENDNA } from './opendna.nf'

params.fasta        = null
params.models       = null
params.lockfile     = null
params.task         = 'disorder'
params.parts        = 10
params.outdir       = 'results'

// 0 means one sequence per embedding batch, which makes the output
// independent of how the input was divided. A larger value is faster on CPU
// and changes the last digit of a few values in ten thousand, so it also makes
// the result depend on the part count. Left at 0 here because a pipeline's
// part count is a scheduling decision and should not be able to change a
// number in a results table.
params.batch_tokens = 0

workflow {
    if (!params.fasta || !params.models || !params.lockfile) {
        error "--fasta, --models and --lockfile are all required"
    }

    OPENDNA(
        file(params.fasta, checkIfExists: true),
        file(params.models, checkIfExists: true),
        file(params.lockfile, checkIfExists: true),
        params.task,
        params.parts
    )
}
