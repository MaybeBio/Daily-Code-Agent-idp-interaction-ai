nextflow.enable.dsl=2

// ==================== Processes ==================== //

include { FASTA_TO_JSON } from './modules/fasta_to_json.nf'
include {SPLIT_JSON} from './modules/split_json.nf'
include { PROTENIX_PREDICT } from './modules/protenix_predict.nf'
include { COLLECT_CHUNKS } from './modules/collect_chunks.nf'
include { PROCESS_MODELS } from './modules/process_models.nf'
include { PLOT_PLDDT } from './modules/plot_plddt.nf'
include { METAPREDICT_DISORDER } from './modules/metapredict_disorder.nf'
include { PLOT_METAPREDICT } from './modules/plot_metapredict.nf'
include { PSAURON_RUN } from './modules/psauron_run.nf'
include { PLOT_PSAURON } from './modules/plot_psauron.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'
include { CHUNK_IPR } from './modules/chunk_ipr.nf'
include { INTERPROSCAN } from './modules/interproscan.nf'
include { CONCAT_IPR } from './modules/concat_ipr.nf'
include { PLOT_INTERPRO } from './modules/plot_interpro.nf'
include { BUILD_SUMMARY_TSV } from './modules/build_summary_tsv.nf'


// ==================== Branch switches ==================== //
// Params given on the command line may arrive as a String ('false') or a
// Boolean (false) depending on the Nextflow version, so always compare
// the lower-cased string form rather than the raw value.

def isOn(value) {
    return value.toString().trim().toLowerCase() in ['true', 'yes', '1']
}

def isOff(value) {
    return value == null || value.toString().trim().toLowerCase() in ['false', 'no', '0', '']
}

// run_interpro is three-way: off, run InterProScan, or plot a supplied XML
def interproMode(value) {
    if ( isOff(value) ) return 'off'
    if ( isOn(value) )  return 'run'
    return 'xml'
}

// Placeholder files stand in for the outputs of disabled branches, so the
// downstream processes keep a fixed input shape. Each has its own name so
// that two placeholders never collide when staged into the same task.
def placeholder(name) {
    return file("${projectDir}/assets/${name}", checkIfExists: true)
}


workflow {

    def run_protenix    = isOn(params.run_protenix)
    def run_metapredict = isOn(params.run_metapredict)
    def run_psauron     = isOn(params.run_psauron)
    def interpro_mode   = interproMode(params.run_interpro)

    if ( !run_protenix && !run_metapredict && !run_psauron && interpro_mode == 'off' ) {
        error "All analysis branches are disabled. Enable at least one of " +
              "--run_protenix, --run_metapredict, --run_psauron or --run_interpro."
    }

    log.info """\
        Branches enabled:
          protenix    : ${run_protenix}
          metapredict : ${run_metapredict}
          psauron     : ${run_psauron}
          interpro    : ${interpro_mode}
        """.stripIndent()

    fasta_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> tuple(fasta, fasta.baseName) }

    // Each plotting branch adds (dataset_name, order, [pngs]) to this list.
    // `order` keeps the report's section order fixed no matter which
    // branches ran or which finished first.
    def report_parts = []

    // ================= METAPREDICT BRANCH (optional) ================= //
    // (run first so its CSV is available for the disorder-coloured
    //  GMM scatter inside PLOT_PLDDT)

    if ( run_metapredict ) {
        metapredict_ch      = METAPREDICT_DISORDER(fasta_ch)
        metapredict_plot_ch = PLOT_METAPREDICT(metapredict_ch)
    } else {
        metapredict_ch      = Channel.empty()
        metapredict_plot_ch = Channel.empty()
    }

    // ================= PROTENIX BRANCH (optional) ================= //

    if ( run_protenix ) {
        json_ch = FASTA_TO_JSON(fasta_ch)

        split_ch  = SPLIT_JSON(json_ch)
        chunks_ch = split_ch.transpose()

        pred_ch = PROTENIX_PREDICT(chunks_ch)

        collected_ch = COLLECT_CHUNKS(pred_ch.groupTuple(by: 0))

        pkl_ch = PROCESS_MODELS(collected_ch)

        // PLOT_PLDDT takes the metapredict CSV for the disorder-coloured
        // scatter. Without metapredict it gets a placeholder and skips
        // just that one plot.
        if ( run_metapredict ) {
            plddt_input_ch = pkl_ch.join(metapredict_ch, by: 0)
        } else {
            def no_metapredict = placeholder('NO_FILE_METAPREDICT')
            plddt_input_ch = pkl_ch.map { dataset_name, pkl -> tuple(dataset_name, pkl, no_metapredict) }
        }
        plddt_plot_ch = PLOT_PLDDT(plddt_input_ch)
    } else {
        chunks_ch     = Channel.empty()
        pkl_ch        = Channel.empty()
        plddt_plot_ch = Channel.empty()
    }

    // ================= PSAURON BRANCH (optional) ================= //

    if ( run_psauron ) {
        psauron_ch      = PSAURON_RUN(fasta_ch)
        psauron_plot_ch = PLOT_PSAURON(psauron_ch)
    } else {
        psauron_ch      = Channel.empty()
        psauron_plot_ch = Channel.empty()
    }

    // ================= INTERPRO BRANCH (optional) ================= //

    if ( interpro_mode == 'xml' ) {

        // User supplied a pre-existing XML — plot only, paired with the
        // dataset's own FASTA so the new script can compute totals.
        def supplied_xml = file(params.run_interpro, checkIfExists: true)
        if ( supplied_xml.extension.toLowerCase() != 'xml' ) {
            error "params.run_interpro must point to an InterProScan XML " +
                  "file (got: ${supplied_xml})"
        }
        interpro_inputs_ch = fasta_ch
            .map { fasta_file, dataset_name ->
                tuple(dataset_name, fasta_file, supplied_xml)
            }
        interpro_plot_ch = PLOT_INTERPRO(interpro_inputs_ch)

    } else if ( interpro_mode == 'run' ) {

        if ( run_protenix ) {
            // Split into the same number of chunks Protenix ended up with,
            // so both branches parallelize the same way.
            ipr_chunk_count_ch = chunks_ch
                .map { dataset_name, chunk_file -> tuple(dataset_name, 1) }
                .groupTuple(by: 0)
                .map { dataset_name, ones -> tuple(dataset_name, ones.size()) }
        } else {
            // No Protenix chunks to mirror, so derive the count from the
            // FASTA using the same chunk_size Protenix would have used.
            ipr_chunk_count_ch = fasta_ch
                .map { fasta_file, dataset_name ->
                    def n_seqs = fasta_file.countFasta()
                    def n_chunks = Math.max(1, (int) Math.ceil(n_seqs / (params.chunk_size as double)))
                    tuple(dataset_name, n_chunks)
                }
        }

        chunk_ipr_input_ch = fasta_ch
            .map { fasta_file, dataset_name -> tuple(dataset_name, fasta_file) }
            .join(ipr_chunk_count_ch, by: 0)
            .map { dataset_name, fasta_file, n_chunks -> tuple(fasta_file, dataset_name, n_chunks) }

        ipr_chunks_ch = CHUNK_IPR(chunk_ipr_input_ch).transpose()

        ipr_results_ch = INTERPROSCAN(ipr_chunks_ch)

        ipr_collected_ch = ipr_results_ch
            .groupTuple(by: 0)
            .map { dataset_name, tsv_list, xml_list ->
                def tsvs = tsv_list instanceof List ? tsv_list.flatten() : [tsv_list]
                def xmls  = xml_list  instanceof List ? xml_list.flatten()  : [xml_list]
                tuple(dataset_name, tsvs, xmls)
            }

        ipr_concat_ch = CONCAT_IPR(ipr_collected_ch)

        // Pair the merged XML back up with each dataset's FASTA for plotting
        interpro_inputs_ch = ipr_concat_ch
            .map { dataset_name, tsv, xml -> tuple(dataset_name, xml) }
            .join(
                fasta_ch.map { fasta_file, dataset_name -> tuple(dataset_name, fasta_file) },
                by: 0
            )
            .map { dataset_name, xml, fasta_file -> tuple(dataset_name, fasta_file, xml) }

        interpro_plot_ch = PLOT_INTERPRO(interpro_inputs_ch)

    } else {

        interpro_plot_ch   = Channel.empty()
        interpro_inputs_ch = Channel.empty()

    }

    // ================= REPORT STAGE ================= //
    // Gather the PNGs from whichever branches ran, keep them in the usual
    // section order (pLDDT, metapredict, psauron, InterPro) and build one
    // PDF per dataset.

    def asList = { x -> x instanceof List ? x : [x] }

    if ( run_protenix ) {
        report_parts << plddt_plot_ch.map { dataset_name, pngs -> tuple(dataset_name, 1, asList(pngs)) }
    }
    if ( run_metapredict ) {
        report_parts << metapredict_plot_ch.map { dataset_name, pngs -> tuple(dataset_name, 2, asList(pngs)) }
    }
    if ( run_psauron ) {
        report_parts << psauron_plot_ch.map { dataset_name, pngs -> tuple(dataset_name, 3, asList(pngs)) }
    }
    if ( interpro_mode != 'off' ) {
        // Both InterPro PNGs go into the report; the TSV is published to
        // results/interpro_plots but not embedded.
        report_parts << interpro_plot_ch.map { dataset_name, bar_png, dist_png, ipr_tsv ->
            tuple(dataset_name, 4, [bar_png, dist_png])
        }
    }

    def n_report_parts = report_parts.size()

    report_ch = Channel.empty()
        .mix(*report_parts)
        // Every enabled branch emits exactly one tuple per dataset, so
        // `size` lets each dataset's report start as soon as it is ready.
        .groupTuple(by: 0, size: n_report_parts)
        .map { dataset_name, orders, png_lists ->
            def pngs = [orders, png_lists].transpose()
                .sort { it[0] }
                .collectMany { it[1] }
            tuple(dataset_name, pngs)
        }

    GENERATE_REPORT(report_ch)

    // ================= SUMMARY TSV STAGE ================= //
    // One-row TSV per dataset with the headline plot metrics. Columns for
    // branches that did not run are left blank so the schema stays fixed.

    summary_ch = fasta_ch.map { fasta_file, dataset_name -> tuple(dataset_name, fasta_file) }

    if ( run_protenix ) {
        summary_ch = summary_ch.join(pkl_ch, by: 0)
    } else {
        def no_plddt = placeholder('NO_FILE_PLDDT')
        summary_ch = summary_ch.map { row -> tuple(*row, no_plddt) }
    }

    if ( run_metapredict ) {
        summary_ch = summary_ch.join(metapredict_ch, by: 0)
    } else {
        def no_metapredict = placeholder('NO_FILE_METAPREDICT')
        summary_ch = summary_ch.map { row -> tuple(*row, no_metapredict) }
    }

    if ( run_psauron ) {
        summary_ch = summary_ch.join(psauron_ch, by: 0)
    } else {
        def no_psauron = placeholder('NO_FILE_PSAURON')
        summary_ch = summary_ch.map { row -> tuple(*row, no_psauron) }
    }

    if ( interpro_mode != 'off' ) {
        summary_ch = summary_ch.join(
            interpro_inputs_ch.map { dataset_name, fasta_file, xml -> tuple(dataset_name, xml) },
            by: 0
        )
    } else {
        def no_ipr = placeholder('NO_FILE_IPR_XML')
        summary_ch = summary_ch.map { row -> tuple(*row, no_ipr) }
    }

    // Resulting tuple: (dataset_name, fasta, pkl, metapredict_csv, psauron_csv, xml)
    BUILD_SUMMARY_TSV(summary_ch)

}
