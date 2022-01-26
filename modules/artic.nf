// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "scheme/${params.schemeDir}/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "scheme/${params.schemeDir}/${params.schemeVersion}/*.scheme.bed" , emit: bed
    path "scheme" , emit: scheme

    script:

    if(params.bed && params.ref){
        bed = file(params.bed)
        ref = file(params.ref)
        """
        mkdir -p scheme/${params.schemeDir}/${params.schemeVersion}/
        cp ${bed} scheme/${params.schemeDir}/${params.schemeVersion}/${params.schemeDir}.scheme.bed
        cp ${reffasta} scheme/${params.schemeDir}/${params.schemeVersion}/${params.schemeDir}.reference.fasta
        """
    } else {
        """
        git clone ${params.schemeRepoURL} scheme
        """
    }
}

process articGuppyPlex {
    tag { params.prefix + "-" + fastqDir }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}*.fastq", mode: "copy"

    input:
    path(fastqDir)

    output:
    path "${params.prefix}*.fastq", emit: fastq

    script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.prefix} \
    --directory ${fastqDir}
    """
}

process articMinIONMedaka {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    input:
    tuple file(fastq), file(schemeRepo)

    output:
    file("${sampleName}.*")
    
    tuple val(sampleName), file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple val(sampleName), file("${sampleName}.sorted.bam"), emit: mapped
    tuple val(sampleName), file("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple val(sampleName), file("${sampleName}.pass.vcf.gz"), emit: vcf

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion --medaka \
    ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo} \
    --read-file ${fastq}\
    --medaka-model ${params.medaka_model}\
    ${params.schemeDir}/${params.schemeVersion} ${sampleName}
    """
}

process articMinIONNanopolish {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    input:
    tuple file(fastq), file(schemeRepo), file(fast5Pass), file(seqSummary)

    output:
    file("${sampleName}.*")
    
    tuple val(sampleName), file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple val(sampleName), file("${sampleName}.sorted.bam"), emit: mapped
    tuple val(sampleName), file("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple val(sampleName), file("${sampleName}.pass.vcf.gz"), emit: vcf

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo} \
    --read-file ${fastq} \
    --fast5-directory ${fast5Pass} \
    --sequencing-summary ${seqSummary} \
    ${params.schemeDir}/${params.schemeVersion} \
    ${sampleName}
    """
}

process articRemoveUnmappedReads {
    tag { sampleName }

    input:
    tuple val(sampleName), path(bamfile)

    output:
    tuple val(sampleName), file("${sampleName}.mapped.sorted.bam")

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.sorted.bam ${bamfile} 
    """
}

