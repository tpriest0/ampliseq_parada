process DADA2_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.svg"), emit: svg
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "versions.yml"               , emit: versions
    path "*.args.txt"                 , emit: args

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "prefix"
    def args = task.ext.args ?: ''
    def seed = task.ext.seed ?: '100'
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE), method = "radix")
        fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE), method = "radix")

        loess_error_function_mod <- function(trans){
            #' This function calculates error estimates for nucleotide transition probabilities using local polynomial regression
            #' fitting (LOESS).
            #' It is a modified version of the dada2 loessErrfun.
            #'
            #' @param trans A transition matrix containing observed transition counts between nucleotides.
            #' @return Returns an error matrix estimating transition probabilities between nucleotides.
            #' The matrix represents the error rates for different nucleotide transitions.
            #'
            #' @details The function iterates through different combinations of nucleotide transitions (e.g., A to C, G to T, etc.)
            #' and computes error estimates using LOESS regression on observed transition counts. It calculates the ratio of
            #' observed errors to total counts and performs LOESS regression on these ratios against the positions in the sequence.
            #' The resulting predictions are used to estimate error rates for different transitions.
            #' Additionally, it applies constraints to ensure error rates fall within defined boundaries.
            
            qq <- as.numeric(colnames(trans))
            est <- matrix(0, nrow = 0, ncol = length(qq))
            for (nti in c("A", "C", "G", "T")) {
                for (ntj in c("A", "C", "G", "T")) {
                    if (nti != ntj) {
                        errs <- trans[paste0(nti, "2", ntj), ]
                        tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
                        rlogp <- log10((errs + 1) / tot)
                        rlogp[is.infinite(rlogp)] <- NA
                        df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
                        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)
                        pred <- predict(mod.lo, qq)
                        maxrli <- max(which(!is.na(pred)))
                	    minrli <- min(which(!is.na(pred)))
                        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
                        pred[seq_along(pred) < minrli] <- pred[[minrli]]
                        est <- rbind(est, 10 ^ pred)
                    }
                }
            }
            MAX_ERROR_RATE <- 0.25
            MIN_ERROR_RATE <- 1e-07
            est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
            est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
            err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ],
                1 - colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12, ]))
            rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
            colnames(err) <- colnames(trans)
            return(err)
        }

        sink(file = "${prefix}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE, errorEstimationFunction = loess_error_function_mod)
        saveRDS(errF, "${prefix}_1.err.rds")
        errR <- learnErrors(fnRs, $args, multithread = $task.cpus, verbose = TRUE, errorEstimationFunction = loess_error_function_mod)
        saveRDS(errR, "${prefix}_2.err.rds")
        sink(file = NULL)

        pdf("${prefix}_1.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}_1.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        pdf("${prefix}_2.err.pdf")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()
        svg("${prefix}_2.err.svg")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}_1.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        sink(file = "${prefix}_2.err.convergence.txt")
        dada2:::checkConvergence(errR)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${prefix}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${prefix}.err.rds")
        sink(file = NULL)

        pdf("${prefix}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}
