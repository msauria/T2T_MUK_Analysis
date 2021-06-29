SYSTEM  = "linux.x86_64"
BINSIZE = 100000
MINSIZE = 5000
GENOMES = ["hg38", "chm13v1"]

GENOME_SET = "|".join(GENOMES)


rule all:
    input:
        expand("plots/{genome}_mu_pairs.png", genome=GENOMES),
        expand("plots/{plot}.pdf", plot=[
            "combined_mu_karyotype", "individual_mu_karyotypes",
            "mu_horizontal_key", "mu_vertical_key"])


######## Get genomic data #########

rule get_chm13v1_fasta:
    output:
        fa="fasta/chm13v1.fasta",
        sizes="fasta/chm13v1.chrom.sizes"
    log:
        "logs/fasta/chm13v1_download.log"
    shell:
        """
        wget -O {output.fa}.gz http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/genome/t2t-chm13-v1.0.fa.gz
        gzip -d -c {output.fa}.gz > {output.fa}
        cat {output.fa} | awk 'BEGIN{{ OFS=""; ORS=""; tab="\t"; nl="\n"; tl=0 }}\
                               {{ if ( NR == 1 ){{ split($1,A,">"); print A[2],tab; }} \
                                  else {{ if ( $1 ~ /^>/ ){{ split($1,A,">"); print tl,nl,A[2],tab; tl = 0;}} \
                                          else {{ tl = tl + length; }} }} }}\
                               END{{ print tl,nl }}' | sort -r -k2,2n > {output.sizes}
        """

rule get_hg38_fasta:
    output:
        fa="fasta/hg38.fasta",
        sizes="fasta/hg38.chrom.sizes"
    log:
        "logs/fasta/hg38_download.log"
    shell:
        """
        wget -O {output.fa}.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
        gzip -d -c {output.fa}.gz > {output.fa}
        cat {output.fa} | awk 'BEGIN{{ OFS=""; ORS=""; tab="\t"; nl="\n"; tl=0 }}\
                               {{ if ( NR == 1 ){{ split($1,A,">"); print A[2],tab; }} \
                                  else {{ if ( $1 ~ /^>/ ){{ split($1,A,">"); print tl,nl,A[2],tab; tl = 0;}} \
                                          else {{ tl = tl + length; }} }} }}\
                               END{{ print tl,nl }}' | sort -r -k2,2n > {output.sizes}
        """

rule get_chm13v1_annotations:
    input:
        "bin/bigBedToBed"
    output:
        "data/chm13v1_cenSat_annotations.bed"
    log:
        "logs/fasta/chm13v1_annotations.log"
    shell:
        """
        wget -O {output}.tmp.bb http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/cenSatAnnotation.bigBed
        {input} {output}.tmp.bb {output}.tmp
        awk 'BEGIN{{ OFS="\t"}}{{ split($4,A,"_"); print $1,$2,$3,A[1] }}' \
            {output}.tmp > {output}
        rm {output}.tmp*
        """

rule get_hg38_annotations:
    output:
        "data/hg38_cenSat_annotations.bed"
    log:
        "logs/fasta/hg38_annotations.log"
    shell:
        """
        wget -O {output}.tmp.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
        gunzip {output}.tmp.gz
        awk 'BEGIN{{A="censat"; OFS="\t"}}{{ print $2,$3,$4,A }}' {output}.tmp > {output}
        rm {output}.tmp*
        """

######## Get software #########

rule download_software:
    output:
        "bin/{sw}"
    wildcard_constraints:
        sw="wigToBigWig|bigBedToBed"
    params:
        sw="{sw}",
        system=SYSTEM
    log:
        "logs/sw/{sw}.log"
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/admin/exe/{params.system}/{params.sw} -O {output}
        chmod a+rx {output}
        """

rule build_minUniqueKmers:
    output:
        "bin/minUniqueKmer/find_minUniqueKmer.sh"
    log:
        "logs/sw/minUniqueKmer.log"
    shell:
        """
        cd bin && git clone https://github.com/msauria/minUniqueKmer.git && cd ..
        chmod a+rx {output}
        """


######## Make MUs #########

rule make_MU_wigs:
    input:
        fa="fasta/{genome}.fasta",
        sw="bin/minUniqueKmer/find_minUniqueKmer.sh"
    output:
        mul="fasta/{genome}.fasta.mul.wig",
        mur="fasta/{genome}.fasta.mur.wig",
        sa="fasta/{genome}.fasta.sa"
    wildcard_constraints:
        genome=GENOME_SET
    log:
        "logs/mu/{genome}_mu.wigs.log"
    shell:
        """
        {input.sw} {input.fa}
        """

rule make_MU_bigwig:
    input:
        mu="fasta/{genome}.fasta.{mu}.wig",
        sizes="fasta/{genome}.chrom.sizes",
        sw="bin/wigToBigWig"
    output:
        "BigWigs/{genome}_{mu}.bw"
    wildcard_constraints:
        genome=GENOME_SET,
        mu="mul|mur"
    log:
        "logs/mu/{genome}_{mu}.bigwigs.log"
    shell:
        """
        {input.sw} {input.mu} {input.sizes} {output}
        """


######## Analysis #########

rule find_mu_hist:
    input:
        "BigWigs/{genome}_mul.bw"
    output:
        "results/{genome}_mu_hist.txt"
    wildcard_constraints:
        genome=GENOME_SET
    log:
        "logs/mu/{genome}.hist.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bin/find_mu_histogram.py {input} {output}
        """

rule find_mu_bg:
    input:
        "BigWigs/{genome}_mul.bw"
    output:
        "results/{genome}_mu_binned.bg"
    params:
        binsize=BINSIZE
    wildcard_constraints:
        genome=GENOME_SET
    log:
        "logs/mu/{genome}.binned.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bin/find_binned_mu.py {input} {output} {params.binsize}
        """

rule find_mu_peaks:
    input:
        "BigWigs/{genome}_mul.bw"
    output:
        "results/{genome}_mu_peaks.txt"
    params:
        minsize=MINSIZE
    wildcard_constraints:
        genome=GENOME_SET
    log:
        "logs/mu/{genome}.peaks.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bin/find_mu_peaks.py {input} {output} {params.minsize}
        """

rule find_mu_seq_pairs:
    input:
        fa="fasta/{genome}.fasta",
        sa="fasta/{genome}.fasta.sa",
        peaks="results/{genome}_mu_peaks.txt"
    output:
        "results/{genome}_mu_pairs.txt"
    params:
        minsize=MINSIZE
    wildcard_constraints:
        genome=GENOME_SET
    log:
        "logs/mu/{genome}.pairs.log"
    shell:
        """
        bin/find_mu_pairs.py {input.fa} {input.sa} {input.peaks} \
            {output} {params.minsize}
        """


######## Plotting #########

rule plot_combined_karyotype:
    input:
        "results/chm13v1_mu_hist.txt",
        "results/hg38_mu_hist.txt",
        "fasta/chm13v1.chrom.sizes",
        "fasta/hg38.chrom.sizes"
    output:
        pdf="plots/combined_mu_karyotype.pdf",
        key="plots/combined_mu_karyotype_key.txt"
    log:
        "logs/R/combined_mu_karyotype.log"
    conda:
        "envs/circos.yaml"
    shell:
        """
        Rscript bin/plot_mul_combined_karyotype.R {input} {output.pdf}
        """

rule plot_karyotype:
    input:
        "results/chm13v1_mu_binned.bg",
        "results/hg38_mu_binned.bg",
        "fasta/chm13v1.chrom.sizes",
        "fasta/hg38.chrom.sizes"
    output:
        pdf="plots/individual_mu_karyotypes.pdf"
    log:
        "logs/R/individual_mu_karyotype.log"
    conda:
        "envs/circos.yaml"
    shell:
        """
        Rscript bin/plot_mul_karyotype.R {input} {output.pdf}
        """

rule plot_key:
    input:
        "plots/combined_mu_karyotype_key.txt"
    output:
        "plots/mu_{orientation}_key.pdf"
    params:
        orientation="{orientation}"
    log:
        "logs/R/{orientation}_key.log"
    conda:
        "envs/circos.yaml"
    shell:
        """
        Rscript bin/plot_key_{params.orientation}.R {input} {output}
        """

rule plot_circos_data:
    input:
        sizes="fasta/{genome}.chrom.sizes",
        anno="data/{genome}_cenSat_annotations.bed",
        pairs="results/{genome}_mu_pairs.txt",
        bg="results/{genome}_mu_binned.bg",
        key="plots/combined_mu_karyotype_key.txt"
    output:
        #karyo="circos_data/{genome}_karyotype.txt",
        #links="circos_data/{genome}_links.txt",
        #missing="circos_data/{genome}_missing_mu.txt",
        #mu="circos_data/{genome}_mu.txt",
        #conf="circos_etc/{genome}_circos.conf",
        plot="plots/{genome}_mu_pairs.png"
    params:
        genome="{genome}",
        minsize=MINSIZE
    log:
        "logs/circos/{genome}.plot.log"
    conda:
        "envs/circos.yaml"
    shell:
        """
        bin/format_circos_data.py {input.sizes} {input.anno} \
            {input.pairs} {input.bg} {input.key} {params.minsize} \
            115000 circos_etc/circos.conf {params.genome}
        circos -conf circos_etc/{params.genome}_circos.conf
        """
