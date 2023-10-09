import "gate_sub.wdl" as sub

task cov_pheno {

    String docker
    File phenolistfile
    File phenomatrix
    File covmatrix

    Int disk_size = ceil(size(phenomatrix, "GB") + size(covmatrix, "GB") + 2)

    command <<<

        set -euxo pipefail
        
        Rscript - <<EOF

        # Combine covariate and phenotype matrices, check that phenos listed exist and there is [pheno]_survTime col also

        suppressPackageStartupMessages(library(data.table))

        cov <- fread("${covmatrix}")
        pheno <- fread("${phenomatrix}")
        phenolist <- readLines("${phenolistfile}")
        phenolist <- c(phenolist, paste0(phenolist, "_survTime"))

        not_found <- setdiff(phenolist, names(pheno))
        if (length(not_found) > 0) stop("Columns not found from phenomatrix: ", paste(not_found, collapse = ", "))

        if (! any(c("FID", "IID", "eid") %in% names(pheno))) {
            stop("No ID column (FID, IID, eid) found from phenomatrix")
        } else if ("eid" %in% names(pheno)) {
            pheno[, FID := eid]
            pheno[, IID := eid]
            pheno[, eid := NULL]
        }

        R <- cov[pheno, , on = c("FID", "IID")]

		fwrite(R, file = "cov_phenomatrix.tsv", sep = "\t", quote = FALSE, na = "NA", col.names = TRUE, row.names = FALSE)

		EOF
        
    >>>

    output {

        File cov_phenomatrix = "cov_phenomatrix.tsv"
    }

    runtime {

        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "us-east1-b us-east1-c us-east1-d"
        noAddress: true
    }
}

task null {

    String pheno
    File bedfile
    File bimfile = sub(bedfile, "\\.bed$", ".bim")
    File famfile = sub(bedfile, "\\.bed$", ".fam")
    String prefix = basename(bedfile, ".bed") + "-" + pheno
    File cov_pheno
    String covariates
    String sampleID
    String traitType
    String docker
    Int cpu
    Int mem
    String loco
    Float eventTimeBinSize
    String pcgforUhatforSurvAnalysis 
    Float traceCVcutoff
    Float ratioCVcutoff
    Int minCovariateCount

    Int disk_size = ceil(size(bedfile, "GB") + size(cov_pheno, "GB") * 3 + 4)

    command {

        step1_fitNULLGLMM.R \
            --plinkFile=${sub(bedfile, "\\.bed$", "")} \
            --phenoFile=${cov_pheno} \
            --phenoCol=${pheno} \
            --eventTimeCol=${pheno}_survTime \
            --eventTimeBinSize=${eventTimeBinSize} \
            --covarColList=${covariates} \
            --sampleIDColinphenoFile=${sampleID} \
            --traitType=${traitType} \
            --outputPrefix=${prefix} \
            --nThreads=${cpu} \
            --LOCO=${loco} \
            --traceCVcutoff=${traceCVcutoff} \
            --ratioCVcutoff=${ratioCVcutoff} \
            --minCovariateCount=${minCovariateCount} \
            --pcgforUhatforSurvAnalysis=${pcgforUhatforSurvAnalysis}
    }

    output {

        File modelfile = prefix + ".rda"
        File varianceratio = prefix + ".varianceRatio.txt"
        File randommarkers = prefix + "_30markers.SAIGE.results.txt"
    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "70 GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "us-east1-b us-east1-c us-east1-d"
        noAddress: true
    }
}

workflow gate {

    String docker
    File phenolistfile
    File phenomatrix
    File covmatrix
    String traitType
    Array[String] phenos = read_lines(phenolistfile)
    String loco

    call cov_pheno {
        input:
            docker = docker,
            phenolistfile = phenolistfile,
            phenomatrix = phenomatrix,
            covmatrix = covmatrix
    }

    scatter (pheno in phenos) {

        call null {
            input:
                docker = docker,
                pheno = pheno,
                cov_pheno = cov_pheno.cov_phenomatrix,
                traitType = traitType,
                loco = loco
        }

        call sub.test_combine {
            input:
                docker = docker,
                pheno = pheno,
                nullfile = null.modelfile
        }
    }
}
