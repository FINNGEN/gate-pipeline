import "gate_sub.wdl" as sub

task null {

    String pheno
    File bedfile
    File bimfile = sub(bedfile, "\\.bed$", ".bim")
    File famfile = sub(bedfile, "\\.bed$", ".fam")
    String prefix = basename(bedfile, ".bed") + "-" + pheno
    File phenofile
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

    command {

        step1_fitNULLGLMM.R \
            --plinkFile=${sub(bedfile, "\\.bed$", "")} \
            --phenoFile=${phenofile} \
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
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        noAddress: true
    }
}

workflow gate {

    String docker
    File phenolistfile
    String traitType
    Array[String] phenos = read_lines(phenolistfile)
    String loco

    scatter (pheno in phenos) {

        call null {
            input: docker=docker, pheno=pheno, traitType=traitType, loco=loco
        }

        call sub.test_combine {
            input: docker=docker, pheno=pheno, nullfile=null.modelfile
        }
    }
}
