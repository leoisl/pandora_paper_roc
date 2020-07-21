params.reference_directory = ""
params.reads_directory = ""
params.sample_id = ""

params.help = false
params.testing = false
params.pipeline_root = ""
params.final_outdir = "."
params.max_forks = 10

if (params.help){
    log.info"""
	Pipeline for running medaka.
        Usage: nextflow run medaka.nf <arguments>
        Required arguments:
          --reference_directory	DIRECTORY	Directory containing different references to use when calling
          --reads_directory	DIRECTORY	Directory containing read fastq[gz] to run on
          --sample_id The sample ID

        Optional:
	  --testing
	  --pipeline_root
	  --final_outdir
	  --max_forks
    """.stripIndent()

    exit 0
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Output directory not found: ${params.final_outdir} -- aborting"
}

if (params.reference_directory) {
    reference_assemblies_in = Channel.fromPath("${params.reference_directory}/*.fna.gz")
}
else {
    exit 1, "Reference assembly directory not provided -- aborting"
}

if (params.reads_directory) {
    reads_directory = file(params.reads_directory).toAbsolutePath()
}
else {
    exit 1, "Expected read files not found"
}

process unzip_reference_assembly {
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    file reference_assembly from reference_assemblies_in

    output:
    file "*.{fa,fna,fasta}" into reference_assemblies_medaka

    """
    v=${reference_assembly}
    if [ \${v: -3} == ".gz" ]
    then
        zcat ${reference_assembly} | awk '{print \$1;}' | tr '[:lower:]' '[:upper:]' | tr URYSWKMBDHV N > \${v::-3}
    else
        cat ${reference_assembly} | awk '{print \$1;}' | tr '[:lower:]' '[:upper:]' | tr URYSWKMBDHV N > n.\$v
	mv n.\$v \$v
    fi
    """
}

Channel.fromPath("${reads_directory}/${params.sample_id}.*.random.nanopore.fastq").set{ nanopore_medaka }

medaka_input = nanopore_medaka.combine(reference_assemblies_medaka)

process medaka {
    memory { 64.GB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
    maxRetries 1
    container {
        'shub://leoisl/Singularity_recipes:medaka'
    }
    input:
    set nanopore_reads, reference_assembly from medaka_input

    publishDir final_outdir, mode: 'copy', overwrite: true

    output:
    set(file("medaka_*.vcf"), file("medaka_*.ref.fa")) into medaka_vcf

    """
    # Following ARTIC pipeline from https://github.com/nanoporetech/medaka/issues/173#issuecomment-642566213 :
    medaka --version
    minimap2 -ax map-ont ${reference_assembly} ${nanopore_reads} | samtools view -b - > nanopore.bam
    samtools sort nanopore.bam > nanopore.sorted.bam
    samtools index nanopore.sorted.bam
    medaka consensus --model r941_min_high_g344 nanopore.sorted.bam medaka_consensus_out.hdf
    medaka variant --verbose ${reference_assembly} medaka_consensus_out.hdf medaka_variants.vcf
    # medaka tools annotate medaka_variants.vcf ${reference_assembly} nanopore.sorted.bam medaka_variants.annotated.vcf # this is commented out because we have a bug (see https://github.com/iqbal-lab/pandora1_paper/issues/121#issuecomment-661879246 )

    v=\$(head -n1 ${reference_assembly})
    ref_id=\${v:1:\${#v}}
    covg=${nanopore_reads}
    covg=\${covg##*/}
    covg=\${covg%.random*}
    covg=\${covg##*.}
    cp medaka_variants.vcf medaka_\$covg.\$ref_id.vcf
    cp ${reference_assembly} medaka_\$covg.\$ref_id.ref.fa
    """
}