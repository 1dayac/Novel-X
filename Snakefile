SAMPLE='CHM13_180GB_CrG_GRCh38_phased_possorted'
GIT_ROOT='~/novel_insertions/'
BLAST_DB='/opt/ncbi-db/data/nt'
READGROUP='40587_MissingLibrary_1_unknown_fc'
GENOME='/home/dmeleshko/genome/genome.fa'

rule all:
    input:
        expand("{sample}.vcf", sample=SAMPLE)

rule extract_unmapped:
    input:
        "sample/{sample}.sorted.bam"
    output:
        "unmapped/{sample}.bam"
    shell:
        """
        #samtools sort -m 30G -n {input} -o sample/{wildcards.sample}.sorted.bam
        ~/bxtools/bin/bxtools filter sample/{wildcards.sample}.sorted.bam -b -s 0.2 -q 10 >{output}
        """

rule convert_bam_to_fastq:
    input:
        "unmapped/{sample}.bam"
    output:
        fastq='reads/{sample}.fastq',
        temp_dir='temp_reads_{sample}'
    shell:
        """
        ~/bamtofastq {input} {output.temp_dir}
        ~/longranger-2.1.6/longranger-cs/2.1.6/bin/longranger basic --id reads --fastqs {output.temp_dir}/*
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
        python ~/novel_insertions/discard_singles.py {input.bam} unmapped/{wildcards.sample}.no_singles.bam
        ~/bamtofastq unmapped/{wildcards.sample}.no_singles.bam temp_reads
        ~/longranger-2.1.6/longranger-cs/2.1.6/bin/longranger basic --id reads_for_velvet_{wildcards.sample} --fastqs temp_reads/*
        gunzip reads_for_velvet_{wildcards.sample}/outs/barcoded.fastq.gz
        bash ~/novel_insertions/interleave.sh < reads_for_velvet_{wildcards.sample}/outs/barcoded.fastq reads_for_velvet_{wildcards.sample}/outs/R1.fastq reads_for_velvet_{wildcards.sample}/outs/R2.fastq
        ~/velvet/velveth velvet_{wildcards.sample} 63 -shortPaired -fastq -separate reads_for_velvet_{wildcards.sample}/outs/R1.fastq reads_for_velvet_{wildcards.sample}/outs/R2.fastq
        ~/velvet/velvetg velvet_{wildcards.sample} -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no
        mkdir fasta
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
        filtered_fasta='filtered/{sample}_filtered.fasta',
        megablast='blast/{sample}.megablast',
        cleanmega='blast/{sample}.cleanmega',
        contaminants='blast/{sample}.contaminants'
    shell:
        """
        blastn -task megablast -query {input.filtered_fasta} -db {BLAST_DB} -num_threads 24 > {output.megablast}
        {GIT_ROOT}/cleanmega {output.megablast} {output.cleanmega}
        {GIT_ROOT}/find_contaminations.py {output.cleanmega} {output.contaminants} 
        python {GIT_ROOT}/remove_contaminations.py {output.contaminants} {input.filtered_fasta} {output.filtered_fasta}
        """



rule align_to_contigs:
    input:
        filtered_fasta='filtered/{sample}_filtered.fasta',
        temp_dir='temp_reads_{sample}'
    output:
        mapped_bam='mapped/{sample}.mapped.bam',
        refdata='refdata-{sample}_filtered'
    shell:
        """
        ~/longranger-2.1.6/longranger-cs/2.1.6/bin/longranger mkref {input.filtered_fasta}
        ~/longranger-2.1.6/longranger-cs/2.1.6/bin/longranger align --id=temp_{wildcards.sample} --reference={output.refdata} --fastqs={input.temp_dir}/{READGROUP}
        samtools view -b -F 12 temp_{wildcards.sample}/outs/possorted_bam.bam >{output.mapped_bam}    
        """

rule extract_barcode_list:
    input:
        mapped_bam='mapped/{sample}.mapped.bam'
    output:
        barcode_folder='{sample}_barcodes'
    shell:
        """
        ~/bxtools/bin/bxtools filter {input.mapped_bam} -s 0.2 -c 0.2 -q 50 >mapped/{wildcards.sample}.filtered.bam
        ~/bxtools/bin/bxtools split-by-ref mapped/{wildcards.sample}.filtered.bam -o {output.barcode_folder}
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
        ~/bxtools/bin/bxtools extract {input.sample} {input.barcode_folder} {output.small_bams}/
        """



rule prepare_reads_for_reassembly:
    input:
         small_bams='small_bams_{sample}'
    output:
         small_reads='small_reads_{sample}',
         temp_small_reads='temp_small_reads_{sample}'
    shell:
         """
         mkdir -p {output.temp_small_reads}
         mkdir -p {output.small_reads}
         function prepare_reads {{
         
         a="$(basename $1 | sed "s/\..*//")"
         if [ -d "{output.small_reads}_2/$a" ]; then
         mv {output.small_reads}_2/$a {output.small_reads}/$a
         return
         fi
         ~/bamtofastq {input.small_bams}/$a.bam {output.temp_small_reads}-$a
         mv {output.temp_small_reads}-$a {output.temp_small_reads}/$a
         ~/longranger-2.1.6/longranger-cs/2.1.6/bin/longranger basic --id {output.small_reads}-$a --fastqs {output.temp_small_reads}/$a/{READGROUP}
         mkdir -p {output.small_reads}/$a
         mv {output.small_reads}-$a/outs/barcoded.fastq.gz {output.small_reads}/$a/
         rm -rf {output.small_reads}-$a
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
        gunzip {input.small_reads}/$a/barcoded.fastq.gz
        bash {GIT_ROOT}/interleave.sh < {input.small_reads}/$a/barcoded.fastq {input.small_reads}/$a/R1.fastq {input.small_reads}/$a/R2.fastq
        python2.7 ~/algorithmic-biology-2/assembler/spades.py --only-assembler -t 1 -k 55 --cov-cutoff 3 --pe1-1 {input.small_reads}/$a/R1.fastq --pe1-2 {input.small_reads}/$a/R2.fastq -o {output.assemblies_folder}/$a
        cp {output.assemblies_folder}/$a/scaffolds.fasta {output.contigs}/$a.fasta
        rm -r {output.assemblies_folder}/$a/K55
        }}
        export -f local_assembly
        parallel --jobs 16 local_assembly ::: {input.small_reads}/*
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
        ~/quast/quast.py --min-contig 199 -R {input.contigs}/$a.fasta {output.splitted_insertions}/$a.fasta -o quast_res/$a
        python {GIT_ROOT}/filter_correct_record.py {input.contigs}/$a.fasta quast_res/$a/contigs_reports/all_alignments_$a.tsv {output.filtered_contigs}/$a.fasta
        }}
        export -f filter_target_contigs
        parallel --jobs 16 filter_target_contigs ::: {input.contigs}/*
        cat {output.filtered_contigs}/*.fasta >{output.contigs}
        rm -r quast_res
        """

rule align_to reference:
    input:
        contigs='final_set_{sample}.fasta',
    output:
        final_alignments='all_alignments_{sample}.tsv'
    shell:
        """
        ~/quast/quast.py -R {GENOME} {input.contigs} -o quast_res
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


