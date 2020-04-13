params.reference_directory = ""
params.reads_directory = ""
params.sample_id = ""
params.coverage = ""

params.help = false
params.testing = false
params.pipeline_root = ""
params.final_outdir = "."
params.max_forks = 10

if (params.help){
    log.info"""
	Pipeline for running pandora and other variant callers and comparing the results.
        Usage: nextflow run compare_genotypers.nf <arguments>
        Required arguments:
          --reference_directory	DIRECTORY	Directory of different references to use when calling with e.g. snippy
	At least one required:
          --reads_directory	DIRECTORY	Directory containing read fastq[gz] to run on
          --sample_id		STRING		Sample ID

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
    reference_assemblies_in = Channel.fromPath("${params.reference_directory}/*.{fa,fa.gz,fasta,fasta.gz,fna,fna.gz}")
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
    file "*.{fa,fna,fasta}" into reference_assemblies_snippy

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

illumina = Channel.fromFilePairs("${reads_directory}/${params.sample_id}.${params.coverage}.random.illumina.{1,2}.fastq", flat:true)
snippy_input = illumina.combine(reference_assemblies_snippy).view()

process snippy_genotype_pe_illumina {
    memory { 24.GB }
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container {
        'shub://rmcolq/Singularity_recipes:snippy'
    }   
    cpus 8
    maxForks 5
    
    publishDir final_outdir, mode: 'copy', overwrite: true
    
    input:
    set name,illumina_reads_1,illumina_reads_2,reference_assembly from snippy_input
    
    output:
    set(file("snippy_*.vcf"), file("snippy_*.ref.fa")) into snippy_vcf
    
    """
    snippy --version
    snippy --cpus 8 --outdir snippy_outdir --reference ${reference_assembly} --pe1 ${illumina_reads_1} --pe2 ${illumina_reads_2}
    v=\$(head -n1 ${reference_assembly})
    ref_id=\${v:1:\${#v}}
    covg=${illumina_reads_1}
    covg=\${covg##*/}
    covg=\${covg%.random*}
    covg=\${covg##*.}
    cp snippy_outdir/snps.filt.vcf snippy_\$covg.\$ref_id.vcf
    cp ${reference_assembly} snippy_\$covg.\$ref_id.ref.fa
    """
}

