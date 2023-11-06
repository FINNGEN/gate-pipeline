version 1.0

task cov_pheno {

    input {
        String docker
        File phenolistfile
        File phenomatrix
        File covmatrix

        Int disk_size = ceil(size(phenomatrix, "GB") + size(covmatrix, "GB") + 2)
    }

    command <<<

        set -euxo pipefail
        
        Rscript - <<EOF

        # Combine covariate and phenotype matrices, check that phenos listed exist and there is [pheno]_survTime col also

        suppressPackageStartupMessages(library(data.table))

        cov <- fread("~{covmatrix}")
        pheno <- fread("~{phenomatrix}")
        phenolist <- readLines("~{phenolistfile}")
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
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk ~{disk_size} HDD"
        zones: "us-east1-b us-east1-c us-east1-d"
        noAddress: true
    }
}

task null {

    input {
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
    }

    command <<<

        set -euxo pipefail

        # Remove NA samples from cov_pheno
        awk 'BEGIN {FS=OFS="\t"; p="~{pheno}"} NR==1 { for(i=1;i<=NF;i++) { h[$i]=i } print } NR>1 && $h[p]!="NA" && $h[p"_survTime"]!="NA" { print }' ~{cov_pheno} > cov_phenomatrix.tsv

        step1_fitNULLGLMM.R \
            --plinkFile=~{sub(bedfile, "\\.bed$", "")} \
            --phenoFile=cov_phenomatrix.tsv \
            --phenoCol=~{pheno} \
            --eventTimeCol=~{pheno}_survTime \
            --eventTimeBinSize=~{eventTimeBinSize} \
            --covarColList=~{covariates} \
            --sampleIDColinphenoFile=~{sampleID} \
            --traitType=~{traitType} \
            --outputPrefix=~{prefix} \
            --nThreads=~{cpu} \
            --LOCO=~{loco} \
            --traceCVcutoff=~{traceCVcutoff} \
            --ratioCVcutoff=~{ratioCVcutoff} \
            --minCovariateCount=~{minCovariateCount} \
            --pcgforUhatforSurvAnalysis=~{pcgforUhatforSurvAnalysis}
    >>>

    output {
        File modelfile = prefix + ".rda"
        File varianceratio = prefix + ".varianceRatio.txt"
        File randommarkers = prefix + "_30markers.SAIGE.results.txt"
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "70 GB"
        disks: "local-disk ~{disk_size} HDD"
        zones: "us-east1-b us-east1-c us-east1-d"
        noAddress: true
    }
}

task test {

    input {
        File nullfile
        File varianceratiofile = sub(nullfile, ".rda", ".varianceRatio.txt")
        String outfileprefix = basename(nullfile, ".rda") + "-"
        Array[File] bgenfiles
        File samplefile
        Int minmac
        String docker
    }

    command <<<

        python3 <<EOF
        import os
        import subprocess
        import time
        processes = set()
        cmd_prefix = 'export MKL_NUM_THREADS=1; export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false; \
            step2_SPAtests.R \
            --minMAC=~{minmac} \
            --sampleFile=~{samplefile} \
            --GMMATmodelFile=~{nullfile} \
            --varianceRatioFile=~{varianceratiofile} \
            --numLinesOutput=1000 \
            --IsOutputAFinCaseCtrl=TRUE '
        for file in '~{sep=" " bgenfiles}'.split(' '):
            cmd = cmd_prefix + '--bgenFile=' + file
            cmd = cmd + ' --SAIGEOutputFile=~{outfileprefix}' + os.path.basename(file) + '.SAIGE.txt'
            logfile = open('SAIGE_log_~{outfileprefix}' + os.path.basename(file) + '.txt', 'w')
            processes.add(subprocess.Popen(cmd, shell=True, stdout=logfile))
        print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(len(processes)) + ' processes started', flush=True)
        n_rc0 = 0
        while n_rc0 < len(processes):
            time.sleep(60)
            n_rc0 = 0
            for p in processes:
                p_poll = p.poll()
                if p_poll is not None and p_poll > 0:
                    raise Exception('subprocess returned ' + str(p_poll))
                if p_poll == 0:
                    n_rc0 = n_rc0 + 1
            print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(n_rc0) + ' processes finished', flush=True)
        EOF
    >>>

    output {
        Array[File] out = glob("*.SAIGE.txt")
        Array[File] logs = glob("SAIGE_log_*.txt")
    }

    runtime {
        docker: "~{docker}"
        cpu: length(bgenfiles)
        memory: (4 * length(bgenfiles)) + " GB"
        disks: "local-disk " + (length(bgenfiles) * ceil(size(bgenfiles[0], "G")) + 5) + " HDD"
        zones: "us-east1-b us-east1-c us-east1-d"
        preemptible: 3
        noAddress: true
    }
}

task combine {

    input {
        String pheno
        Array[Array[File]] results2D
        Array[File] results = flatten(results2D)
        String chrcol
        String af_col
        String p_valcol
        String bp_col
        Int loglog_pval
        String docker
        String prefix
        Boolean logP
        String logPStr = if logP then "True" else "False"
    }

    command <<<

        set -e

        echo "`date` concatenating results to ~{prefix}~{pheno}.saige.gz"
        cat \
        <(head -n 1 ~{results[0]}) \
        <(for file in ~{sep=" " results}; do tail -n+2 $file; done) \
        | bgzip > ~{prefix}~{pheno}.saige.gz

        echo "`date` converting results to ~{prefix}~{pheno}.gz"

        python3 <<EOF | sort -k 1,1g -k 2,2g | bgzip > ~{prefix}~{pheno}.gz

        import math, gzip
        from collections import OrderedDict
        from functools import reduce

        def conv_base(logp):
            return logp / math.log(10)

        def red(obj, func_list):
            return "NA" if obj == "NA" else reduce(lambda o, func: func[0](o, *func[1]) if func[0] is not str.format else func[0](func[1], o), func_list, obj)

        mapping = OrderedDict([
            ("#chrom", ("CHR", [(str.replace, ("chr", "")), (str.replace, ("X", "23")), (str.replace, ("Y", "24")), (str.replace, ("MT", "25")), (str.replace, ("M", "25"))])),
            ("pos", ("POS", [(float, ()), (int, ())])), # convert 1e6 to 1000000
            ("ref", ("Allele1", [])),
            ("alt", ("Allele2", [])),
            ("pval", ("p.value", [(float, ()), (math.exp, ()), (str.format, ("{:.2e}"))])) if ~{logPStr} else ("pval", ("p.value", [(float, ()), (str.format, ("{:.2e}"))])),
            ("mlogp", ("p.value", [(float, ()), (abs, ()), (conv_base, ()), (str.format, ("{:.4f}"))])) if ~{logPStr} else ("mlogp", ("p.value", [(float, ()), (math.log10, ()), (abs, ()), (str.format, ("{:.4f}"))])),
            ("beta", ("BETA", [(float, ()), (str.format, ("{:.5f}"))])),
            ("sebeta", ("SE", [(float, ()), (str.format, ("{:.5f}"))])),
            ("af_alt", ("AF_Allele2", [(float, ()), (str.format, ("{:.2e}"))])),
            ("af_alt_cases", ("AF.Events", [(float, ()), (str.format, ("{:.2e}"))])),
            ("af_alt_controls", ("AF.Censored", [(float, ()), (str.format, ("{:.2e}"))])),
        ])

        with gzip.open("~{prefix}~{pheno}.saige.gz", 'rt') as f:
            header = {h:i for i,h in enumerate(f.readline().strip().split(' '))}
            for col in [v[0] for v in mapping.values()]:
                if col not in header.keys():
                    raise Exception('column ' + col + ' not in given file')
            print('\t'.join(mapping.keys()))
            for line in f:
                s = line.strip().split(' ')
                print('\t'.join(str(red(s[header[v[0]]], v[1])) for v in mapping.values()))
        
        EOF

        echo "`date` Manhattan plot start"
        /plot_scripts/ManhattanPlot.r --input=~{prefix}~{pheno}.gz  --PVAL="~{p_valcol}" --knownRegionFlank=1000000 --prefix="~{prefix}~{pheno}"  --ismanhattanplot=TRUE --isannovar=FALSE --isqqplot=FALSE --CHR="~{chrcol}" --POS="~{bp_col}" --ALLELE1=ref --ALLELE2=alt

        echo "`date` QQ plot start"
        /plot_scripts/QQplot.r --input=~{prefix}~{pheno}.gz  --prefix="~{prefix}~{pheno}" --af="~{af_col}" --pvalue="~{p_valcol}"

        echo "`date` tabixing"
        tabix -S 1 -b 2 -e 2 -s 1 ~{prefix}~{pheno}.gz
        echo "`date` done"
        
    >>>

    output {
        File saige_out = prefix + pheno + ".saige.gz"
        File out = prefix + pheno + ".gz"
        File out_ind = prefix + pheno + ".gz.tbi"
        File out_regions = prefix + pheno + ".regions.txt"
        File out_tophits = prefix + pheno + ".tophits.txt" 
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "~{docker}"
        cpu: 1
        memory: "20 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-b us-east1-c us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow gate {

    input {
        String docker
        File phenolistfile
        File phenomatrix
        File covmatrix
        String traitType
        Array[String] phenos = read_lines(phenolistfile)
        String loco
        File bgenlistfile
        Array[Array[String]] bgenfiles2D = read_tsv(bgenlistfile)
    }

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

        scatter (bgenfiles in bgenfiles2D) {
            call test {
                input:
                    docker = docker,
                    nullfile = null.modelfile,
                    bgenfiles = bgenfiles
            }
        }

        call combine {
            input:
                pheno = pheno,
                results2D = test.out
        }
    }

    output {
            Array[File] nullmodels = null.modelfile
            Array[File] sumstats = combine.out
            Array[File] sumstats_tbi = combine.out
            Array[File] regions = combine.out_regions
            Array[File] tophits = combine.out_tophits
            Array[Array[File]] pngs = combine.pngs
    }
}
