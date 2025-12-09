nextflow.enable.dsl=2

params.fastq1 = params.fastq1 ?: null
params.fastq2 = params.fastq2 ?: null
params.reference = params.reference
params.prefix = params.prefix ?: 'sample'
params.outdir = params.outdir ?: "${baseDir}/results"
params.simulate_reads = params.simulate_reads?.toBoolean() ?: false


// workflow
workflow {

    //force reference
    if( !params.reference ) {
        error "You must provide --reference"
    }

    //if simulated reads
    if( params.simulate_reads ) {
        Channel
            .fromPath(params.reference)
            .set { reference_fasta }
        sim_output = mutate_and_simulate(reference_fasta)
        sim_output
            .map { mut_csv, fq1, fq2 -> tuple(fq1, fq2) }
            .set { fastq_pair }
    }

    //if not simulated reads
    else {
        if( !params.fastq1 || !params.fastq2 ) {
            error "simulate_reads = false: you must provide --fastq1 and --fastq2"
        }
        
        fastq_pair = Channel.of([ params.fastq1, params.fastq2 ])}

    alignment_results = alignment(fastq_pair, params.reference)
    variants = variant_calling(alignment_results, params.reference)
    comparison = compare_vcfs(variants)
}


// processes
process mutate_and_simulate {
    executor 'local'
    tag { infile }

    input:
    path infile

    output:
    tuple path("${params.prefix}_mutations.csv"),
          path("${params.prefix}_fastq1.fq"),
          path("${params.prefix}_fastq2.fq")

    """
    mkdir -p ${params.outdir}

    python ${projectDir}/scripts/simulate_reads.py \
        -r ${infile} \
        -n ${params.prefix} \
        -o .

    cp ${params.prefix}_mutations.csv ${params.outdir}/
    cp ${params.prefix}_fastq1.fq      ${params.outdir}/
    cp ${params.prefix}_fastq2.fq      ${params.outdir}/
    """
}

process alignment {
    executor 'local'
    tag { params.prefix }

    input:
    tuple path(fq1), path(fq2)
    val  reference_fasta

    output:
    tuple path("${params.prefix}.bam"),
          path("${params.prefix}.bam.bai")

    """
    mkdir -p ${params.outdir}

    minimap2 -a -x sr ${projectDir}/${reference_fasta} ${fq1} ${fq2} \
        | samtools view -h -F 0x900 - \
        | samtools sort -O bam \
        > ${params.prefix}.bam

    samtools index ${params.prefix}.bam

    cp ${params.prefix}.bam     ${params.outdir}/
    cp ${params.prefix}.bam.bai ${params.outdir}/
    """
}

process variant_calling {
    executor 'local'
    tag { params.prefix }

    input:
    tuple path(bam), path(bai)
    val   reference_fasta

    output:
    tuple path("${params.prefix}_bcftools.vcf"),
          path("${params.prefix}_snippy.vcf")

    """
    mkdir -p ${params.outdir}

    bcftools mpileup -Ou \
      -f ${projectDir}/${reference_fasta} \
      -Q 0 -A -d 1000000 ${bam} | \
    bcftools call -mv -Ov -o ${params.prefix}_bcftools.vcf

    snippy \
      --outdir unknown \
      --ref ${projectDir}/${reference_fasta} \
      --bam ${bam} \
      --prefix variants \
      --mapqual 0 \
      --mincov 1 \
      --minfrac 0.1


    cp unknown/variants.vcf ${params.prefix}_snippy.vcf
    cp ${params.prefix}_bcftools.vcf ${params.outdir}/
    cp ${params.prefix}_snippy.vcf   ${params.outdir}/
    """
}

process compare_vcfs {
    executor 'local'
    tag { params.prefix }

    input:
    tuple path(bcf_vcf), path(snippy_vcf)

    output:
    tuple path("${params.prefix}_compare.txt"),
          path("${params.prefix}_compare_summary.txt")

    """
    mkdir -p ${params.outdir}

    TRUTH_FILE="${params.outdir}/${params.prefix}_mutations.csv"

    python ${projectDir}/scripts/compare_vcfs.py \
        -b ${bcf_vcf} \
        -s ${snippy_vcf} \
        -t ${params.outdir}/${params.prefix}_mutations.csv \
        -o ${params.prefix}

    cp ${params.prefix}_compare.txt         ${params.outdir}/
    cp ${params.prefix}_compare_summary.txt ${params.outdir}/
    """
}
