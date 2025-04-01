/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT AND VARIABLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Input

if (params.metadata) {
    ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
} else { ch_metadata = Channel.empty() }

if (params.dada_ref_tax_custom) {
    //custom ref taxonomy input from params.dada_ref_tax_custom & params.dada_ref_tax_custom_sp
    ch_assigntax = Channel.fromPath("${params.dada_ref_tax_custom}", checkIfExists: true)
    if (params.dada_ref_tax_custom_sp) {
        ch_addspecies = Channel.fromPath("${params.dada_ref_tax_custom_sp}", checkIfExists: true)
    } else { ch_addspecies = Channel.empty() }
    ch_dada_ref_taxonomy = Channel.empty()
    val_dada_ref_taxonomy = "user"
} else if (params.dada_ref_taxonomy && !params.skip_dada_taxonomy && !params.skip_taxonomy) {
    //standard ref taxonomy input from params.dada_ref_taxonomy & conf/ref_databases.config
    ch_dada_ref_taxonomy = Channel.fromList(params.dada_ref_databases[params.dada_ref_taxonomy]["file"]).map { file(it) }
    val_dada_ref_taxonomy = params.dada_ref_taxonomy.replace('=','_').replace('.','_')
} else {
    ch_dada_ref_taxonomy = Channel.empty()
    val_dada_ref_taxonomy = "none"
}

if (params.sintax_ref_taxonomy && !params.skip_taxonomy) {
    ch_sintax_ref_taxonomy = Channel.fromList(params.sintax_ref_databases[params.sintax_ref_taxonomy]["file"]).map { file(it) }
    val_sintax_ref_taxonomy = params.sintax_ref_taxonomy.replace('=','_').replace('.','_')
} else {
    ch_sintax_ref_taxonomy = Channel.empty()
    val_sintax_ref_taxonomy = "none"
}

// report sources
ch_report_template = Channel.fromPath("${params.report_template}", checkIfExists: true)
ch_report_css = Channel.fromPath("${params.report_css}", checkIfExists: true)
ch_report_logo = Channel.fromPath("${params.report_logo}", checkIfExists: true)
ch_report_abstract = params.report_abstract ? Channel.fromPath(params.report_abstract, checkIfExists: true) : []

// Set non-params Variables

single_end = params.single_end
if (params.pacbio || params.iontorrent) {
    single_end = true
}

trunclenf = params.trunclenf ?: 0
trunclenr = params.trunclenr ?: 0
if ( !single_end && !params.illumina_pe_its && (params.trunclenf == null || params.trunclenr == null) && !params.input_fasta ) {
    find_truncation_values = true
    log.warn "No DADA2 cutoffs were specified (`--trunclenf` & `--trunclenr`), therefore reads will be truncated where median quality drops below ${params.trunc_qmin} (defined by `--trunc_qmin`) but at least a fraction of ${params.trunc_rmin} (defined by `--trunc_rmin`) of the reads will be retained.\nThe chosen cutoffs do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\n"
} else { find_truncation_values = false }

// save params to values to be able to overwrite it
tax_agglom_min = params.tax_agglom_min
tax_agglom_max = params.tax_agglom_max

//use custom taxlevels from --dada_assign_taxlevels or database specific taxlevels if specified in conf/ref_databases.config
if ( params.dada_ref_taxonomy ) {
    taxlevels = params.dada_assign_taxlevels ? "${params.dada_assign_taxlevels}" :
        params.dada_ref_databases[params.dada_ref_taxonomy]["taxlevels"] ?: ""
} else { taxlevels = params.dada_assign_taxlevels ? "${params.dada_assign_taxlevels}" : "" }
if ( params.sintax_ref_taxonomy ) {
    sintax_taxlevels = params.sintax_ref_databases[params.sintax_ref_taxonomy]["taxlevels"] ?: ""
} else {
    sintax_taxlevels = ""
}

//make sure that taxlevels adheres to requirements when mixed with addSpecies
if ( params.dada_ref_taxonomy && !params.skip_dada_addspecies && !params.skip_dada_taxonomy && !params.skip_taxonomy && taxlevels ) {
    if ( !taxlevels.endsWith(",Genus,Species") && !taxlevels.endsWith(",Genus") ) {
        error("Incompatible settings: To use exact species annotations, taxonomic levels must end with `,Genus,Species` or `,Genus` but are currently `${taxlevels}`. Taxonomic levels can be set with `--dada_assign_taxlevels`. Skip exact species annotations with `--skip_dada_addspecies`.\n")
    }
}

// This tracks tax tables produced during pipeline and each table will be used during phyloseq
ch_tax_for_phyloseq = Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE & SUBWORKFLOW: Installed directly from nf-core/modules & nf-core/subworkflows
//

include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { VSEARCH_CLUSTER                   } from '../modules/nf-core/vsearch/cluster/main'
include { FASTA_NEWICK_EPANG_GAPPA          } from '../subworkflows/nf-core/fasta_newick_epang_gappa/main'

//
// MODULE: Installed directly from nf-core/modules
//

include { RENAME_RAW_DATA_FILES         } from '../modules/local/rename_raw_data_files'
include { DADA2_ERR                     } from '../modules/local/dada2_err'
include { NOVASEQ_ERR                   } from '../modules/local/novaseq_err'
include { DADA2_DENOISING               } from '../modules/local/dada2_denoising'
include { DADA2_RMCHIMERA               } from '../modules/local/dada2_rmchimera'
include { DADA2_STATS                   } from '../modules/local/dada2_stats'
include { DADA2_MERGE                   } from '../modules/local/dada2_merge'
include { DADA2_SPLITREGIONS            } from '../modules/local/dada2_splitregions'
include { BARRNAP                       } from '../modules/local/barrnap'
include { BARRNAPSUMMARY                } from '../modules/local/barrnapsummary'
include { FILTER_SSU                    } from '../modules/local/filter_ssu'
include { FILTER_LEN as FILTER_LEN_ASV  } from '../modules/local/filter_len'
include { FILTER_LEN as FILTER_LEN_ITSX } from '../modules/local/filter_len'
include { MERGE_STATS as MERGE_STATS_FILTERSSU    } from '../modules/local/merge_stats'
include { MERGE_STATS as MERGE_STATS_FILTERLENASV } from '../modules/local/merge_stats'
include { MERGE_STATS as MERGE_STATS_CODONS       } from '../modules/local/merge_stats'
include { FILTER_CODONS                 } from '../modules/local/filter_codons'
include { FORMAT_FASTAINPUT             } from '../modules/local/format_fastainput'
include { FORMAT_TAXONOMY               } from '../modules/local/format_taxonomy'
include { ITSX_CUTASV                   } from '../modules/local/itsx_cutasv'
include { MERGE_STATS as MERGE_STATS_STD} from '../modules/local/merge_stats'
include { FORMAT_PPLACETAX              } from '../modules/local/format_pplacetax'
include { FILTER_STATS                  } from '../modules/local/filter_stats'
include { MERGE_STATS as MERGE_STATS_FILTERTAXA } from '../modules/local/merge_stats'
include { METADATA_ALL                  } from '../modules/local/metadata_all'
include { METADATA_PAIRWISE             } from '../modules/local/metadata_pairwise'
include { SBDIEXPORT                    } from '../modules/local/sbdiexport'
include { SBDIEXPORTREANNOTATE          } from '../modules/local/sbdiexportreannotate'
include { SUMMARY_REPORT                } from '../modules/local/summary_report'
include { FILTER_CLUSTERS               } from '../modules/local/filter_clusters'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PARSE_INPUT                   } from '../subworkflows/local/parse_input'
include { DADA2_PREPROCESSING           } from '../subworkflows/local/dada2_preprocessing'
include { CUTADAPT_WORKFLOW             } from '../subworkflows/local/cutadapt_workflow'
include { DADA2_TAXONOMY_WF             } from '../subworkflows/local/dada2_taxonomy_wf'
include { SINTAX_TAXONOMY_WF            } from '../subworkflows/local/sintax_taxonomy_wf'

//
// FUNCTIONS
//
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ampliseq_pipeline'
include { makeComplement         } from '../subworkflows/local/utils_nfcore_ampliseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMPLISEQ {

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Create input channels
    //
    ch_input_fasta = Channel.empty()
    ch_input_reads = Channel.empty()
    if ( params.input ) {
        // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
        ch_input_reads = Channel.fromSamplesheet("input")
            .map{ meta, readfw, readrv ->
                meta.single_end = single_end.toBoolean()
                def reads = single_end ? readfw : [readfw,readrv]
                if ( !meta.single_end && !readrv ) { error("Entry `reverseReads` is missing in $params.input for $meta.id, either correct the samplesheet or use `--single_end`, `--pacbio`, or `--iontorrent`") } // make sure that reverse reads are present when single_end isnt specified
                if ( !meta.single_end && ( readfw.getSimpleName() == meta.id || readrv.getSimpleName() == meta.id ) ) { error("Entry `sampleID` cannot be identical to simple name of `forwardReads` or `reverseReads`, please change `sampleID` in $params.input for sample $meta.id") } // sample name and any file name without extensions arent identical, because rename_raw_data_files.nf would forward 3 files (2 renamed +1 input) instead of 2 in that case
                if ( meta.single_end && ( readfw.getSimpleName() == meta.id+"_1" || readfw.getSimpleName() == meta.id+"_2" ) ) { error("Entry `sampleID`+ `_1` or `_2` cannot be identical to simple name of `forwardReads`, please change `sampleID` in $params.input for sample $meta.id") } // sample name and file name without extensions arent identical, because rename_raw_data_files.nf would forward 2 files (1 renamed +1 input) instead of 1 in that case
                return [meta, reads] }
    } else if ( params.input_fasta ) {
        ch_input_fasta = Channel.fromPath(params.input_fasta, checkIfExists: true)
    } else if ( params.input_folder ) {
        PARSE_INPUT ( params.input_folder, single_end, params.multiple_sequencing_runs, params.extension )
        ch_input_reads = PARSE_INPUT.out.reads
    } else {
        error("One of `--input`, `--input_fasta`, `--input_folder` must be provided!")
    }

    //
    // Add primer info to sequencing files
    //
    if ( params.multiregion ) {
        // is multiple region analysis
        ch_input_reads
            .combine( Channel.fromSamplesheet("multiregion") )
            .map{ info, reads, multi ->
                def meta = info + multi
                return [ meta, reads ] }
            .map{ info, reads ->
                def meta = info +
                    [id: info.sample+"_"+info.fw_primer+"_"+info.rv_primer] +
                    [fw_primer_revcomp: makeComplement(info.fw_primer.reverse())] +
                    [rv_primer_revcomp: makeComplement(info.rv_primer.reverse())]
                return [ meta, reads ] }
            .set { ch_input_reads }
    } else {
        // is single region
        ch_input_reads
            .map{ info, reads ->
                def meta = info +
                    [region: null, region_length: null] +
                    [fw_primer: params.FW_primer, rv_primer: params.RV_primer] +
                    [id: info.sample] +
                    [fw_primer_revcomp: params.FW_primer ? makeComplement(params.FW_primer.reverse()) : null] +
                    [rv_primer_revcomp: params.RV_primer ? makeComplement(params.RV_primer.reverse()) : null]
                return [ meta, reads ] }
            .set { ch_input_reads }
    }

    //Filter empty files
    ch_input_reads.dump(tag:'ch_input_reads')
        .branch {
            failed: it[0].single_end ? it[1].countFastq() < params.min_read_counts : it[1][0].countFastq() < params.min_read_counts || it[1][1].countFastq() < params.min_read_counts
            passed: true
        }
        .set { ch_reads_result }
    ch_reads_result.passed.set { ch_reads }
    ch_reads_result.failed
        .map { meta, reads -> [ meta.id ] }
        .collect()
        .subscribe {
            samples = it.join("\n")
            if (params.ignore_empty_input_files) {
                log.warn "At least one input file for the following sample(s) had too few reads (<$params.min_read_counts):\n$samples\nThe threshold can be adjusted with `--min_read_counts`. Ignoring failed samples and continue!\n"
            } else {
                error("At least one input file for the following sample(s) had too few reads (<$params.min_read_counts):\n$samples\nEither remove those samples, adjust the threshold with `--min_read_counts`, or ignore that samples using `--ignore_empty_input_files`.")
            }
        }
    ch_reads.dump(tag: 'ch_reads')

    //
    // MODULE: Rename files
    //
    RENAME_RAW_DATA_FILES ( ch_reads )
    ch_versions = ch_versions.mix(RENAME_RAW_DATA_FILES.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC ( RENAME_RAW_DATA_FILES.out.fastq )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //
    // MODULE: Cutadapt
    //
    if (!params.skip_cutadapt) {
        CUTADAPT_WORKFLOW (
            RENAME_RAW_DATA_FILES.out.fastq,
            params.illumina_pe_its,
            params.double_primer
        ).reads.set { ch_trimmed_reads }
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_WORKFLOW.out.logs.collect{it[1]})
        ch_versions = ch_versions.mix(CUTADAPT_WORKFLOW.out.versions)
    } else {
        ch_trimmed_reads = RENAME_RAW_DATA_FILES.out.fastq
    }

    //
    // SUBWORKFLOW: Read preprocessing & QC plotting with DADA2
    //
    DADA2_PREPROCESSING (
        ch_trimmed_reads,
        single_end,
        find_truncation_values,
        trunclenf,
        trunclenr
    ).reads.set { ch_filt_reads }
    ch_versions = ch_versions.mix(DADA2_PREPROCESSING.out.versions)

    //
    // MODULES: ASV generation with DADA2
    //

    //run error model
    if ( !params.illumina_novaseq ) {
        DADA2_ERR ( ch_filt_reads )
        ch_errormodel = DADA2_ERR.out.errormodel
        ch_versions = ch_versions.mix(DADA2_ERR.out.versions)
    } else {
        DADA2_ERR ( ch_filt_reads )
        NOVASEQ_ERR ( DADA2_ERR.out.errormodel )
        ch_errormodel = NOVASEQ_ERR.out.errormodel
        ch_versions = ch_versions.mix(DADA2_ERR.out.versions)
    }

    //group by meta
    ch_filt_reads
        .join( ch_errormodel )
        .set { ch_derep_errormodel }
    DADA2_DENOISING ( ch_derep_errormodel.dump(tag: 'into_denoising')  )
    ch_versions = ch_versions.mix(DADA2_DENOISING.out.versions)

    DADA2_RMCHIMERA ( DADA2_DENOISING.out.seqtab )
    ch_versions = ch_versions.mix(DADA2_RMCHIMERA.out.versions)

    //group by sequencing run & group by meta
    DADA2_PREPROCESSING.out.logs
        .join( DADA2_DENOISING.out.denoised )
        .join( DADA2_DENOISING.out.mergers )
        .join( DADA2_RMCHIMERA.out.rds )
        .set { ch_track_numbers }
    DADA2_STATS ( ch_track_numbers )
    ch_versions = ch_versions.mix(DADA2_STATS.out.versions)

    //merge if several runs, otherwise just publish
    DADA2_MERGE (
        DADA2_STATS.out.stats.map { meta, stats -> stats }.collect(),
        DADA2_RMCHIMERA.out.rds.map { meta, rds -> rds }.collect() )
    ch_versions = ch_versions.mix(DADA2_MERGE.out.versions)

    //merge cutadapt_summary and dada_stats files
    if (!params.skip_cutadapt) {
        MERGE_STATS_STD (CUTADAPT_WORKFLOW.out.summary, DADA2_MERGE.out.dada2stats)
        ch_stats = MERGE_STATS_STD.out.tsv
        ch_versions = ch_versions.mix(MERGE_STATS_STD.out.versions)
    } else {
        ch_stats = DADA2_MERGE.out.dada2stats
    }

    ch_dada2_fasta = DADA2_MERGE.out.fasta
    ch_dada2_asv = DADA2_MERGE.out.asv

    // MODULE : ASV post-clustering with VSEARCH
    //
    if (params.vsearch_cluster && !params.multiregion) {
        ch_fasta_for_clustering = ch_dada2_fasta
            .map {
                fasta ->
                    def meta = [:]
                    meta.id = "ASV_post_clustering"
                    [ meta, fasta ] }
        VSEARCH_CLUSTER ( ch_fasta_for_clustering )
        ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions)
        FILTER_CLUSTERS ( VSEARCH_CLUSTER.out.clusters, ch_dada2_asv )
        ch_versions = ch_versions.mix(FILTER_CLUSTERS.out.versions)
        ch_dada2_fasta = FILTER_CLUSTERS.out.fasta
        ch_dada2_asv = FILTER_CLUSTERS.out.asv
    }

    //
    // Entry for ASV fasta files via "--input_fasta"
    //
    if ( params.input_fasta ) {
        FORMAT_FASTAINPUT( ch_input_fasta )
        ch_unfiltered_fasta = FORMAT_FASTAINPUT.out.fasta
        ch_versions = ch_versions.mix(FORMAT_FASTAINPUT.out.versions)
    } else {
        ch_unfiltered_fasta = ch_dada2_fasta
    }

    //
    // Modules : Filter rRNA
    //
    if ( !params.skip_barrnap && params.filter_ssu && !params.multiregion ) {
        BARRNAP ( ch_unfiltered_fasta )
        ch_versions = ch_versions.mix(BARRNAP.out.versions)
        BARRNAPSUMMARY ( BARRNAP.out.gff.collect() )
        ch_versions = ch_versions.mix(BARRNAPSUMMARY.out.versions)
        BARRNAPSUMMARY.out.warning.subscribe {
            if ( it.baseName.toString().startsWith("WARNING") ) {
                error("Barrnap could not identify any rRNA in the ASV sequences! This will result in all sequences being removed with SSU filtering.")
            }
        }
        ch_barrnapsummary = BARRNAPSUMMARY.out.summary
        FILTER_SSU ( ch_unfiltered_fasta, ch_dada2_asv.ifEmpty( [] ), BARRNAPSUMMARY.out.summary )
        ch_versions = ch_versions.mix(FILTER_SSU.out.versions)
        MERGE_STATS_FILTERSSU ( ch_stats, FILTER_SSU.out.stats )
        ch_versions = ch_versions.mix(MERGE_STATS_FILTERSSU.out.versions)
        ch_stats = MERGE_STATS_FILTERSSU.out.tsv
        ch_dada2_fasta = FILTER_SSU.out.fasta
        ch_dada2_asv = FILTER_SSU.out.asv
    } else if ( !params.skip_barrnap && !params.filter_ssu && !params.multiregion ) {
        BARRNAP ( ch_unfiltered_fasta )
        BARRNAPSUMMARY ( BARRNAP.out.gff.collect() )
        BARRNAPSUMMARY.out.warning.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn "Barrnap could not identify any rRNA in the ASV sequences. We recommended to use the --skip_barrnap option for these sequences." }
        ch_barrnapsummary = BARRNAPSUMMARY.out.summary
        ch_versions = ch_versions.mix(BARRNAP.out.versions)
        ch_versions = ch_versions.mix(BARRNAPSUMMARY.out.versions)
        ch_dada2_fasta = ch_unfiltered_fasta
    } else {
        ch_barrnapsummary = Channel.empty()
        ch_dada2_fasta = ch_unfiltered_fasta
    }

    //
    // Modules : amplicon length filtering
    //
    if ( (params.min_len_asv || params.max_len_asv) && !params.multiregion ) {
        FILTER_LEN_ASV ( ch_dada2_fasta, ch_dada2_asv.ifEmpty( [] ) )
        ch_versions = ch_versions.mix(FILTER_LEN_ASV.out.versions)
        MERGE_STATS_FILTERLENASV ( ch_stats, FILTER_LEN_ASV.out.stats )
        ch_stats = MERGE_STATS_FILTERLENASV.out.tsv
        ch_dada2_fasta = FILTER_LEN_ASV.out.fasta
        ch_dada2_asv = FILTER_LEN_ASV.out.asv
        // Make sure that not all sequences were removed
        ch_dada2_fasta.subscribe { if (it.countLines() == 0) error("ASV length filtering activated by '--min_len_asv' or '--max_len_asv' removed all ASVs, please adjust settings.") }
    }

    //
    // Modules : Filtering based on codons in an open reading frame
    //
    if ( params.filter_codons && !params.multiregion ) {
        FILTER_CODONS ( ch_dada2_fasta, ch_dada2_asv.ifEmpty( [] ) )
        ch_versions = ch_versions.mix(FILTER_CODONS.out.versions)
        MERGE_STATS_CODONS( ch_stats, FILTER_CODONS.out.stats )
        ch_versions = ch_versions.mix(MERGE_STATS_CODONS.out.versions)
        ch_stats = MERGE_STATS_CODONS.out.tsv
        ch_dada2_fasta = FILTER_CODONS.out.fasta
        ch_dada2_asv = FILTER_CODONS.out.asv
        // Make sure that not all sequences were removed
        ch_dada2_fasta.subscribe { if (it.countLines() == 0) error("ASV codon filtering activated by '--filter_codons' removed all ASVs, please adjust settings.") }
    }

    //
    // Modules : ITSx - cut out ITS region if long ITS reads
    //
    ch_full_fasta = ch_dada2_fasta
    if (params.cut_its == "none") {
        ch_fasta = ch_dada2_fasta
    } else {
        if (params.cut_its == "full") {
            outfile = params.its_partial ? "ASV_ITS_seqs.full_and_partial.fasta" : "ASV_ITS_seqs.full.fasta"
        }
        else if (params.cut_its == "its1") {
            outfile =  params.its_partial ? "ASV_ITS_seqs.ITS1.full_and_partial.fasta" : "ASV_ITS_seqs.ITS1.fasta"
        }
        else if (params.cut_its == "its2") {
            outfile =  params.its_partial ? "ASV_ITS_seqs.ITS2.full_and_partial.fasta" : "ASV_ITS_seqs.ITS2.fasta"
        }
        ITSX_CUTASV ( ch_full_fasta, outfile )
        ch_versions = ch_versions.mix(ITSX_CUTASV.out.versions)
        FILTER_LEN_ITSX ( ITSX_CUTASV.out.fasta, [] )
        ch_fasta = FILTER_LEN_ITSX.out.fasta
    }

    //
    // SUBWORKFLOW / MODULES : Taxonomic classification with DADA2 or SINTAX
    //

    //DADA2
    if (!params.skip_taxonomy && !params.skip_dada_taxonomy) {
        if (!params.dada_ref_tax_custom) {
            //standard ref taxonomy input from conf/ref_databases.config
            FORMAT_TAXONOMY ( ch_dada_ref_taxonomy.collect(), val_dada_ref_taxonomy )
            ch_versions = ch_versions.mix(FORMAT_TAXONOMY.out.versions)
            ch_assigntax = FORMAT_TAXONOMY.out.assigntax
            ch_addspecies = FORMAT_TAXONOMY.out.addspecies
        }
        DADA2_TAXONOMY_WF (
            ch_assigntax,
            ch_addspecies,
            val_dada_ref_taxonomy,
            ch_fasta,
            ch_full_fasta,
            taxlevels
        ).tax.set { ch_dada2_tax }
        ch_versions = ch_versions.mix(DADA2_TAXONOMY_WF.out.versions)
        ch_tax_for_phyloseq = ch_tax_for_phyloseq.mix ( ch_dada2_tax.map { it = [ "dada2", file(it) ] } )
    } else {
        ch_dada2_tax = Channel.empty()
    }

    // SINTAX
    if (!params.skip_taxonomy && params.sintax_ref_taxonomy) {
        SINTAX_TAXONOMY_WF (
            ch_sintax_ref_taxonomy.collect(),
            val_sintax_ref_taxonomy,
            ch_fasta,
            ch_full_fasta,
            sintax_taxlevels
        ).tax.set { ch_sintax_tax }
        ch_versions = ch_versions.mix(SINTAX_TAXONOMY_WF.out.versions)
        ch_tax_for_phyloseq = ch_tax_for_phyloseq.mix ( ch_sintax_tax.map { it = [ "sintax", file(it) ] } )
    } else {
        ch_sintax_tax = Channel.empty()
    }

    // Phylo placement
    if ( params.pplace_tree ) {
        ch_pp_data = ch_fasta.map { it ->
            [ meta: [ id: params.pplace_name ?: 'user_tree' ],
            data: [
                alignmethod:  params.pplace_alnmethod ?: 'hmmer',
                queryseqfile: it,
                refseqfile:   file( params.pplace_aln, checkIfExists: true ),
                hmmfile:      [],
                refphylogeny: file( params.pplace_tree, checkIfExists: true ),
                model:        params.pplace_model,
                taxonomy:     params.pplace_taxonomy ? file( params.pplace_taxonomy, checkIfExists: true ) : []
            ] ]
        }
        FASTA_NEWICK_EPANG_GAPPA ( ch_pp_data )
        ch_versions = ch_versions.mix( FASTA_NEWICK_EPANG_GAPPA.out.versions )
        ch_pplace_tax = FORMAT_PPLACETAX ( FASTA_NEWICK_EPANG_GAPPA.out.taxonomy_per_query ).tsv
        ch_tax_for_phyloseq = ch_tax_for_phyloseq.mix ( PHYLOSEQ_INTAX_PPLACE ( ch_pplace_tax ).tsv.map { it = [ "pplace", file(it) ] } )
    } else {
        ch_pplace_tax = Channel.empty()
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )

        ch_multiqc_report_list = MULTIQC.out.report.toList()
    } else {
        ch_multiqc_report_list = Channel.empty()
    }

    //
    // MODULE: Summary Report
    //
    if (!params.skip_report && !params.multiregion) {
        SUMMARY_REPORT (
            ch_report_template,
            ch_report_css,
            ch_report_logo,
            ch_report_abstract,
            ch_metadata.ifEmpty( [] ),
            params.input ? file(params.input) : [], // samplesheet input
            ch_input_fasta.ifEmpty( [] ), // fasta input
            !params.input_fasta && !params.skip_fastqc && !params.skip_multiqc ? MULTIQC.out.plots : [], //.collect().flatten().collectFile(name: "fastqc_per_sequence_quality_scores_plot.svg")
            !params.skip_cutadapt ? CUTADAPT_WORKFLOW.out.summary.collect().ifEmpty( [] ) : [],
            find_truncation_values,
            DADA2_PREPROCESSING.out.args.first().ifEmpty( [] ),
            !params.skip_dada_quality ? DADA2_PREPROCESSING.out.qc_svg.ifEmpty( [] ) : [],
            !params.skip_dada_quality ? DADA2_PREPROCESSING.out.qc_svg_preprocessed.ifEmpty( [] ) : [],
            DADA2_ERR.out.svg
                .map {
                    meta_old, svgs ->
                    def meta = [:]
                    meta.single_end = meta_old.single_end
                    [ meta, svgs, meta_old.run ] }
                .groupTuple(by: 0 )
                .map {
                    meta_old, svgs, runs ->
                    def meta = [:]
                    meta.single_end = meta_old.single_end
                    meta.run = runs.flatten()
                    [ meta, svgs.flatten() ]
                }.ifEmpty( [[],[]] ),
            DADA2_MERGE.out.asv.ifEmpty( [] ),
            ch_unfiltered_fasta.ifEmpty( [] ), // this is identical to DADA2_MERGE.out.fasta if !params.input_fasta
            DADA2_MERGE.out.dada2asv.ifEmpty( [] ),
            DADA2_MERGE.out.dada2stats.ifEmpty( [] ),
            params.vsearch_cluster ? FILTER_CLUSTERS.out.asv.ifEmpty( [] ) : [],
            !params.skip_barrnap ? BARRNAPSUMMARY.out.summary.ifEmpty( [] ) : [],
            params.filter_ssu ? FILTER_SSU.out.stats.ifEmpty( [] ) : [],
            params.filter_ssu ? FILTER_SSU.out.fasta.ifEmpty( [] ) : [],
            params.min_len_asv || params.max_len_asv ? FILTER_LEN_ASV.out.stats.ifEmpty( [] ) : [],
            params.min_len_asv || params.max_len_asv ? FILTER_LEN_ASV.out.len_orig.ifEmpty( [] ) : [],
            params.filter_codons ? FILTER_CODONS.out.fasta.ifEmpty( [] ) : [],
            params.filter_codons ? FILTER_CODONS.out.stats.ifEmpty( [] ) : [],
            params.cut_its != "none" ? ITSX_CUTASV.out.summary.ifEmpty( [] ) : [],
            !params.skip_taxonomy && params.dada_ref_taxonomy && !params.skip_dada_taxonomy ? ch_dada2_tax.ifEmpty( [] ) : [],
            !params.skip_taxonomy && params.dada_ref_taxonomy && !params.skip_dada_taxonomy ? DADA2_TAXONOMY_WF.out.cut_tax.ifEmpty( [[],[]] ) : [[],[]],
            !params.skip_taxonomy && params.sintax_ref_taxonomy ? ch_sintax_tax.ifEmpty( [] ) : [],
            !params.skip_taxonomy && params.pplace_tree ? ch_pplace_tax.ifEmpty( [] ) : [],
            !params.skip_taxonomy && params.pplace_tree ? FASTA_NEWICK_EPANG_GAPPA.out.heattree.ifEmpty( [[],[]] ) : [[],[]],
        )
        ch_versions    = ch_versions.mix(SUMMARY_REPORT.out.versions)
    }

    //
    // Save input files in results folder
    //
    if ( params.input ) {
        file("${params.outdir}/input").mkdir()
        file("${params.input}").copyTo("${params.outdir}/input")
    }
    if ( params.input_fasta ) {
        file("${params.outdir}/input").mkdir()
        file("${params.input_fasta}").copyTo("${params.outdir}/input")
    }
    if ( params.multiregion ) {
        file("${params.outdir}/input").mkdir()
        file("${params.multiregion}").copyTo("${params.outdir}/input")
    }
    if ( params.metadata ) {
        file("${params.outdir}/input").mkdir()
        file("${params.metadata}").copyTo("${params.outdir}/input")
    }

    emit:
    multiqc_report = ch_multiqc_report_list      // MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
