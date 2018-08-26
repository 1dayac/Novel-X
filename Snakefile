
configfile: "config.json"
configfile: "path_to_executables_config.json"

SAMPLE=config["sample"]
GIT_ROOT=config["root"]
BLAST_DB=config["blast_db"]
READGROUP=config["readgroup"]
GENOME=config["genome"]
SAMTOOLS=config["samtools"]
VELVETH=config["velveth"]
VELVETG=config["velvetg"]
BLASTN=config["blastn"]
BLASTDB=config["blast_db"]
QUAST=config["quast"]
SPADES=config["spades.py"]
LONGRANGER=config["longranger"]
THREADS=cinfig["threads"]
MEMORY=cinfig["memory"]
MEMORY_PER_THREAD=cinfig["memory_per_thread"]
ADDITIONAL_BAMTOFASTQ_FLAGS=cinfig["additional_flags"]
rule all:
    input:
        expand("{sample}.vcf", sample=SAMPLE)

rule extract_unmapped:
    input:
        "sample/{sample}.bam"
    output:
        unmapped="unmapped/{sample}.bam",
        sorted="sample/{sample}.sorted.bam"
    shell:
        """
        {SAMTOOLS} sort -@ {THREADS} -m {MEMORY_PER_THREAD}G -n {input} -o sample/{wildcards.sample}.sorted.bam
        {GIT_ROOT}/bxtools/bin/bxtools filter sample/{wildcards.sample}.sorted.bam -b -s 0.2 -q 10 >{output}
        """

rule convert_bam_to_fastq:
    input:
        "unmapped/{sample}.bam"
    output:
        fastq='reads/{sample}.fastq',
        temp_dir='temp_reads_{sample}'
    shell:
        """
        {GIT_ROOT}/bamtofastq {ADDITIONAL_BAMTOFASTQ_FLAGS} {input} {output.temp_dir}
        {LONGRANGER} basic --id reads --fastqs {output.temp_dir}/*
        mv reads/outs/barcoded.fastq.gz {output.fastq}.gz
        gunzip {output.fastq}.gz
        """


rule deinterleave:
    input:
        "reads/{sample}.fastq"
    output:
        left='reads/R1_{sample}.fastq',
        right='reads/R2_{sample}.fastq'
    shell:
        """
        bash {GIT_ROOT}/interleave.sh < {input} {output.left} {output.right}
        rm {input}
        """

rule velvet_assembly:
    input:
        bam='unmapped/{sample}.bam'
    output:
        fasta='fasta/{sample}.fasta',
        no_singles='unmapped/{sample}.no_singles.bam'
    shell:
        """
        rm -rf temp_reads
        python {GIT_ROOT}/discard_singles.py {input.bam} unmapped/{wildcards.sample}.no_singles.bam
        mkdir temp_reads
        {GIT_ROOT}/bxtools/bin/bxtools bamtofastq unmapped/{wildcards.sample}.no_singles.bam temp_reads/
        {VELVETH} velvet_{wildcards.sample} 63 -shortPaired -fastq -separate temp_reads/{wildcards.sample}.no_singles_R1.fastq temp_reads/{wildcards.sample}.no_singles_R2.fastq
        {VELVETG} velvet_{wildcards.sample} -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no
        rm -rf temp_reads/
        mkdir -p fasta
        cp velvet_{wildcards.sample}/contigs.fa fasta/{wildcards.sample}.fasta
        """


rule filter_length:
    input:
         fasta='fasta/{sample}.fasta'
    output:
         filtered_fasta='filtered/{sample}_filtered.long.fasta'
    shell:
         """
         {GIT_ROOT}/contig_length_filter.py 200 {input.fasta} {output.filtered_fasta}
         """

rule filter_contaminants:
    input:
        filtered_fasta='filtered/{sample}_filtered.long.fasta'
    output:
        filtered_fasta='filtered/{sample}_filtered.fasta'
    run:
        if BLAST_DB != 'None':
            shell("{BLASTN} -task megablast -query {input.filtered_fasta} -db {BLAST_DB} -num_threads {THREADS} > blast/{wildcards.sample}.megablast")
            shell("{GIT_ROOT}/cleanmega blast/{wildcards.sample}.megablast blast/{wildcards.sample}.cleanmega")
            shell("{GIT_ROOT}/find_contaminations.py blast/{wildcards.sample}.cleanmega blast/{wildcards.sample}.contaminants")
            shell("python {GIT_ROOT}/remove_contaminations.py blast/{wildcards.sample}.contaminants {input.filtered_fasta} {output.filtered_fasta}")"
        else:
            shell("cp {input.filtered_fasta} {output.filtered_fasta}")"



rule align_to_contigs:
    input:
        filtered_fasta='filtered/{sample}_filtered.fasta',
        temp_dir='temp_reads_{sample}'
    output:
        mapped_bam='mapped/{sample}.mapped.bam',
        refdata='refdata-{sample}_filtered'
    shell:
        """
        {LONGRANGER} mkref {input.filtered_fasta}
        {LONGRANGER} align --localcores={THREADS} --localmem={MEMORY} --id=temp_{wildcards.sample} --reference={output.refdata} --fastqs={input.temp_dir}/{READGROUP}
        {SAMTOOLS} view -b -F 12 temp_{wildcards.sample}/outs/possorted_bam.bam >{output.mapped_bam}
        """

rule extract_barcode_list:
    input:
        mapped_bam='mapped/{sample}.mapped.bam'
    output:
        barcode_folder='{sample}_barcodes'
    shell:
        """
        {GIT_ROOT}/bxtools/bin/bxtools filter {input.mapped_bam} -s 0.2 -c 0.2 -q 50 >mapped/{wildcards.sample}.filtered.bam
        {GIT_ROOT}/bxtools/bin/bxtools split-by-ref mapped/{wildcards.sample}.filtered.bam -o {output.barcode_folder}
        """

rule assemble_barcode_list:
    input:
        barcode_folder='{sample}_barcodes',
        sample='sample/{sample}.sorted.bam'
    output:
        small_bams='small_bams_{sample}'
    shell:
        """
        mkdir small_bams_{wildcards.sample}
        {GIT_ROOT}/bxtools/bin/bxtools extract {input.sample} {input.barcode_folder} {output.small_bams}/
        """



rule prepare_reads_for_reassembly:
    input:
         small_bams='small_bams_{sample}'
    output:
         small_reads='small_reads_{sample}'
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
        assemblies_folder='local_assemblies_{sample}',
        contigs='contigs_{sample}'
    shell:
        """
        mkdir -p {output.contigs}
        mkdir -p {output.assemblies_folder}
        function local_assembly {{
        a="$(basename $1 | sed "s/\..*q//")"
        python2.7 {SPADES} --only-assembler -t 1 -m {MEMORY_PER_THREAD} -k 77 --cov-cutoff 3 --pe1-1 {input.small_reads}/$a/$a_R1.fastq --pe1-2 {input.small_reads}/$a/$a_R2.fastq -o {output.assemblies_folder}/$a
        cp {output.assemblies_folder}/$a/scaffolds.fasta {output.contigs}/$a.fasta
        rm -rf {output.assemblies_folder}/$a/K55
        }}
        export -f local_assembly
        parallel --jobs {THREADS} local_assembly ::: {input.small_reads}/*
        """


rule filter_target_contigs:
    input:
        contigs='contigs_{sample}',
        insertions='filtered/{sample}_filtered.fasta'
    output:
        filtered_contigs='filtered_contigs_{sample}',
        splitted_insertions='splitted_insertions_{sample}',
        contigs='final_set_{sample}.fasta'
    shell:
        """
        mkdir -p {output.splitted_insertions}
        python {GIT_ROOT}/split_fasta.py {input.insertions} {output.splitted_insertions}
        mkdir -p quast_res
        mkdir -p {output.filtered_contigs}
        function filter_target_contigs {{
        a="$(basename $1 | sed "s/\..*//")"
        {QUAST} -t 1 --min-contig 199 -R {input.contigs}/$a.fasta {output.splitted_insertions}/$a.fasta -o quast_res/$a
        python {GIT_ROOT}/filter_correct_record.py {input.contigs}/$a.fasta quast_res/$a/contigs_reports/all_alignments_$a.tsv {output.filtered_contigs}/$a.fasta
        }}
        export -f filter_target_contigs
        parallel --jobs {THREADS} filter_target_contigs ::: {input.contigs}/*
        cat {output.filtered_contigs}/*.fasta >{output.contigs}
        rm -rf quast_res
        """

rule align_to reference:
    input:
        contigs='final_set_{sample}.fasta',
    output:
        final_alignments='all_alignments_{sample}.tsv'
    shell:
        """
        {QUAST} -R {GENOME} {input.contigs} -o quast_res
        cp quast_res/contigs_reports/all_alignments_final_set_{wildcards.sample}.tsv {output.final_alignments}
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


