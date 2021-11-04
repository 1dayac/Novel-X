
configfile: "config.json"
configfile: "path_to_executables_config.json"

SAMPLE=config["sample"]
GIT_ROOT=config["root"]
BLAST_DB=config["blast_db"]
GENOME=config["genome"]
SAMTOOLS=config["samtools"]
VELVETH=config["velveth"]
VELVETG=config["velvetg"]
BLASTN=config["blastn"]
BLASTDB=config["blast_db"]
QUAST=config["quast"]
SPADES=config["spades.py"]
LONGRANGER=config["longranger"]
THREADS=config["threads"]
MEMORY=config["memory"]
MEMORY_PER_THREAD=config["memory_per_thread"]
ADDITIONAL_BAMTOFASTQ_FLAGS=config["additional_flags"]
SPADES_K=config["spades_k_assembly"]
VELVET_K=config["velvet_k_assembly"]
VELVET_COVERAGE=config["velvet_coverage"]
DATA=config["tenx"]
rule all:
    input:
        expand("{sample}.vcf", sample=SAMPLE)

rule extract_unmapped_reads:
    input:
        "sample/{sample}.bam"
    output:
        unmapped_bam="unmapped_bam/{sample}.bam",
        sorted="sample/{sample}.sorted.bam"
    shell:
        """
        {SAMTOOLS} sort -@ {THREADS} -n {input} -o {output.sorted}
        {GIT_ROOT}/bxtools/bin/bxtools filter {output.sorted} -b -s 0.2 -q 10 >{output.unmapped_bam}
        """

rule assemble_unmapped_reads:
    input:
        bam='unmapped_bam/{sample}.bam'
    output:
        fasta='fasta/{sample}.fasta'
    run:
        shell("rm -rf temp_reads")
        shell("mkdir temp_reads")
        shell("mkdir -p fasta")

        shell("{GIT_ROOT}/bxtools/bin/bxtools bamtofastq {input.bam} temp_reads/")
        if DATA == "other":
            shell("{SPADES} -t {THREADS} -k {VELVET_K} -1 temp_reads/{wildcards.sample}_R1.fastq -2 temp_reads/{wildcards.sample}_R2.fastq --only-assembler -o spades_{wildcards.sample}")
            shell("cp spades_{wildcards.sample}/scaffolds.fasta fasta/{wildcards.sample}.fasta")
        else:
            shell("{VELVETH} velvet_{wildcards.sample} {VELVET_K} -shortPaired -fastq -separate temp_reads/{wildcards.sample}_R1.fastq temp_reads/{wildcards.sample}_R2.fastq")
            shell("{VELVETG} velvet_{wildcards.sample} -exp_cov auto -cov_cutoff {VELVET_COVERAGE} -max_coverage 100 -scaffolding no")
            shell("cp velvet_{wildcards.sample}/contigs.fa fasta/{wildcards.sample}.fasta")
        shell("rm -rf temp_reads/")


rule filter_short_contigs:
    input:
         fasta='fasta/{sample}.fasta'
    output:
         filtered_fasta='filtered/{sample}_filtered.long.fasta'
    shell:
         """
         {GIT_ROOT}/contig_length_filter.py 200 {input.fasta} {output.filtered_fasta}
         """

rule convert_unmapped_bam_to_fastq:
    input:
        "unmapped_bam/{sample}.bam"
    output:
        temp_dir=directory('temp_reads_{sample}')
    shell:
        """
        {GIT_ROOT}/bamtofastq {ADDITIONAL_BAMTOFASTQ_FLAGS} {input} {output.temp_dir}
        """

rule filter_contaminants:
    input:
        filtered_fasta='filtered/{sample}_filtered.long.fasta'
    output:
        filtered_fasta='filtered/{sample}_filtered.fasta'
    run:
        if BLAST_DB != 'None':
            shell("mkdir -p blast")
            shell("{BLASTN} -task megablast -query {input.filtered_fasta} -db {BLAST_DB} -num_threads {THREADS} > blast/{wildcards.sample}.megablast")
            shell("{GIT_ROOT}/cleanmega blast/{wildcards.sample}.megablast blast/{wildcards.sample}.cleanmega")
            shell("{GIT_ROOT}/find_contaminations.py blast/{wildcards.sample}.cleanmega blast/{wildcards.sample}.contaminants")
            shell("python {GIT_ROOT}/remove_contaminations.py blast/{wildcards.sample}.contaminants {input.filtered_fasta} {output.filtered_fasta}")
        else:
            shell("cp {input.filtered_fasta} {output.filtered_fasta}")



rule align_unmapped_reads_to_the_contigs:
    input:
        filtered_fasta='filtered/{sample}_filtered.fasta',
        temp_dir='temp_reads_{sample}'
    output:
        mapped_bam='mapped/{sample}.mapped.bam',
        refdata=directory('refdata-{sample}_filtered')
    shell:
        """
        {LONGRANGER} mkref {input.filtered_fasta}
        {LONGRANGER} align --localcores={THREADS} --localmem={MEMORY} --id=temp_{wildcards.sample} --reference={output.refdata} --fastqs={input.temp_dir}/
        {SAMTOOLS} view -b -F 12 temp_{wildcards.sample}/outs/possorted_bam.bam >{output.mapped_bam}
        """

rule extract_barcode_lists:
    input:
        mapped_bam='mapped/{sample}.mapped.bam'
    output:
        barcode_folder=directory('{sample}_barcodes')
    shell:
        """
        {GIT_ROOT}/bxtools/bin/bxtools filter {input.mapped_bam} -s 0.2 -c 0.2 -q 50 >mapped/{wildcards.sample}.filtered.bam
        {GIT_ROOT}/bxtools/bin/bxtools split-by-ref mapped/{wildcards.sample}.filtered.bam -o {output.barcode_folder}
        """

rule extract_bam_subsets:
    input:
        barcode_folder='{sample}_barcodes',
        sample='sample/{sample}.sorted.bam'
    output:
        small_bams=directory('small_bams_{sample}')
    shell:
        """
        mkdir small_bams_{wildcards.sample}
        {GIT_ROOT}/bxtools/bin/bxtools extract {input.sample} {input.barcode_folder} {output.small_bams}/
        """



rule prepare_reads_for_local_assembly:
    input:
         small_bams='small_bams_{sample}'
    output:
         small_reads=directory('small_reads_{sample}')
    shell:
         """
         mkdir -p {output.small_reads}
         function prepare_reads {{
         
         a="$(basename $1 | sed "s/\..*//")"
         if [ -d "{output.small_reads}_2/$a" ]; then
         mv {output.small_reads}_2/$a {output.small_reads}/$a
         return
         fi
         mkdir -p {output.small_reads}/$a
         {GIT_ROOT}/bxtools/bin/bxtools bamtofastq {input.small_bams}/$a.bam {output.small_reads}/$a/
         }}
         export -f prepare_reads
         parallel --jobs 16 prepare_reads ::: {input.small_bams}/*
         """

rule local_assembly:
    input:
        small_reads='small_reads_{sample}'
    output:
        assemblies_folder=directory('local_assemblies_{sample}'),
        contigs=directory('contigs_{sample}')
    shell:
        """
        mkdir -p {output.contigs}
        mkdir -p {output.assemblies_folder}
        function local_assembly {{
        a="$(basename $1 | sed "s/\..*q//")"
        if [ -f "{output.assemblies_folder}/$a/scaffolds.fasta" ]; then
        cp {output.assemblies_folder}/$a/scaffolds.fasta {output.contigs}/$a.fasta
        return
        fi
        {SPADES} --only-assembler -t 1 -m {MEMORY_PER_THREAD} -k {SPADES_K} --cov-cutoff 3 --pe1-1 {input.small_reads}/$a/$a\_R1.fastq --pe1-2 {input.small_reads}/$a/$a\_R2.fastq -o {output.assemblies_folder}/$a
        cp {output.assemblies_folder}/$a/scaffolds.fasta {output.contigs}/$a.fasta
        rm -rf {output.assemblies_folder}/$a/K*
        }}
        export -f local_assembly
        parallel --jobs {THREADS} local_assembly ::: {input.small_reads}/*
        """


rule filter_target_contigs:
    input:
        contigs='contigs_{sample}',
        insertions='filtered/{sample}_filtered.fasta'
    output:
        filtered_contigs=directory('filtered_contigs_{sample}'),
        prefilter_contigs=directory('prefiltered_contigs_{sample}'),
        splitted_insertions=directory('splitted_insertions_{sample}'),
        contigs='final_set_{sample}.fasta'
    shell:
        """
        mkdir -p {output.splitted_insertions}
        mkdir -p {output.prefilter_contigs}
        python {GIT_ROOT}/split_fasta.py {input.insertions} {output.splitted_insertions}
        mkdir -p quast_res
        mkdir -p {output.filtered_contigs}
        function filter_target_contigs {{
        a="$(basename $1 | sed "s/\..*//")"
        {GIT_ROOT}/contig_length_filter.py 500 {input.contigs}/$a.fasta {output.prefilter_contigs}/$a.fasta
        {QUAST} -t 1 --fast --min-contig 199 -R {output.prefilter_contigs}/$a.fasta {output.splitted_insertions}/$a.fasta -o quast_res/$a
        python {GIT_ROOT}/filter_correct_record.py {output.prefilter_contigs}/$a.fasta quast_res/$a/contigs_reports/all_alignments_$a.tsv {output.filtered_contigs}/$a.fasta
        }}
        export -f filter_target_contigs
        set +oe pipefail
        parallel --jobs {THREADS} filter_target_contigs ::: {input.contigs}/*
        cat {output.filtered_contigs}/*.fasta >{output.contigs}
        rm -rf quast_res
        rm -f core*
        """

rule align_to_reference:
    input:
        contigs='final_set_{sample}.fasta',
    output:
        final_alignments='all_alignments_{sample}.tsv'
    shell:
        """
        {QUAST} --fast -R {GENOME} {input.contigs} -o quast_res
        cp quast_res/contigs_reports/all_alignments_final_set_*.tsv {output.final_alignments}
        """

rule produce_vcf:
    input:
        final_alignments='all_alignments_{sample}.tsv',
        contigs='final_set_{sample}.fasta'
    output:
        final_vcf='{sample}.vcf'
    shell:
        """
        python {GIT_ROOT}/minimap_to_vcf.py {input.final_alignments} {input.contigs} {output.final_vcf}
        """


