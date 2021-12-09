process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz") optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    input:
        path(ref)

    output:
        tuple path('ref.fa'), path('ref.mmi')

    script:
        """
        ln -s ${ref} ref.fa
        minimap2 -x sr -d ref.mmi ref.fa
        """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using minimap2 (https://github.com/lh3/minimap2)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: 'copy'

    input:
        tuple val(sampleName), path(forward), path(reverse), path(ref), path("*")

    output:
        tuple val(sampleName), path("${sampleName}.sorted.bam")

    script:
        """
        minimap2 -a ${ref} ${forward} ${reverse} | \
            samtools view -S -b | \
            samtools sort -o ${sampleName}.sorted.bam
        """
}

process lenFilter {

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.lenfiltered.bam", mode: 'copy'

    input:
        tuple val(sampleName), path("${sampleName}.sorted.bam")
    
    output:
        tuple val(sampleName), path("${sampleName}.lenfiltered.bam")

    script:
        """
        samtools view -h ${sampleName}.sorted.bam | \
            awk 'length(\$10) >= ${params.MinReadLen} && length(\$10) <= ${params.MaxReadLen} || \$1 ~ /^@/' | \
            samtools view -b -o ${sampleName}.lenfiltered.bam
        """     
}



process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bedfile)

    output:
    tuple val(sampleName), path("${sampleName}.mapped.bam"), emit: mapped
    tuple val(sampleName), path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
   
    if ( params.cleanBamHeader )
        """
        samtools reheader --no-PG  -c 'sed "s/${sampleName}/sample/g"' ${bam} | \
        samtools view -F4 -o sample.mapped.bam

        mv sample.mapped.bam ${sampleName}.mapped.bam
        
        samtools index ${sampleName}.mapped.bam

        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out

        samtools reheader --no-PG  -c 'sed "s/${sampleName}/sample/g"' ivar.out.bam | \
        samtools sort -o sample.mapped.primertrimmed.sorted.bam

        mv sample.mapped.primertrimmed.sorted.bam ${sampleName}.mapped.primertrimmed.sorted.bam
        """

    else
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
}

// process callVariants {

//     tag { sampleName }

//     publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

//     input:
//     tuple val(sampleName), path(bam), path(ref)

//     output:
//     tuple val(sampleName), path("${sampleName}.variants.tsv"), emit: variants

//     script:
//         """
//         samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
//         ivar variants -r ${ref} -m ${params.ivarMinDepth} -p ${sampleName}.variants -q ${params.ivarMinVariantQuality} -t ${params.ivarMinFreqThreshold}
//         """
// }

process callVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.norm.vcf", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(ref)

    output:
    val(sampleName), emit: sampleName
    path(ref), emit: ref
    path("${sampleName}.ambiguous.norm.vcf.gz*"), emit: amb_vcf_gz
    path("${sampleName}.mask.txt"), emit: mask
    path("${sampleName}.fixed.norm.vcf.gz*"), emit: fixed_vcf_gz

    script:
        """


        freebayes -p 1 \
                  -f ${ref} \
                  -F 0.2 \
                  -C 1 \
                  --pooled-continuous \
                  --min-coverage ${params.MinDepth} \
                  --gvcf --gvcf-dont-use-chunk true ${bam} > ${sampleName}.gvcf

        # make depth mask, split variants into ambiguous/consensus
        # NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
        process_gvcf.py -d ${params.MinDepth} \
                        -l ${params.MinFreqThreshold} \
                        -u ${params.FreqThreshold} \
                        -m ${sampleName}.mask.txt \
                        -v ${sampleName}.variants.vcf \
                        -c ${sampleName}.consensus.vcf ${sampleName}.gvcf

        # normalize variant records into canonical VCF representation
        for v in "variants" "consensus"; do
            bcftools norm -f ${ref} ${sampleName}.\$v.vcf > ${sampleName}.\$v.norm.vcf
        done

        # split the consensus sites file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF
        for vt in "ambiguous" "fixed"; do
            cat ${sampleName}.consensus.norm.vcf | awk -v vartag=ConsensusTag=\$vt '\$0 ~ /^#/ || \$0 ~ vartag' > ${sampleName}.\$vt.norm.vcf
            bgzip -f ${sampleName}.\$vt.norm.vcf
            tabix -f -p vcf ${sampleName}.\$vt.norm.vcf.gz
        done
        """

}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        val(sampleName)
        path(ref)
        path(amb_vcf_gz)
        path(mask)
        path(fixed_vcf_gz)

    output:
        tuple val(sampleName), path("${sampleName}.primertrimmed.consensus.fa")

    script:
        """
        # apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
        bcftools consensus -f ${ref} -I ${sampleName}.ambiguous.norm.vcf.gz > ${sampleName}.ambiguous.fasta
        # apply remaning variants, including indels
        bcftools consensus -f ${sampleName}.ambiguous.fasta -m ${mask} ${sampleName}.fixed.norm.vcf.gz | sed s/MN908947.3/${sampleName}/ > ${sampleName}.primertrimmed.consensus.fa
        """
}

process cramToFastq {
    /**
    * Converts CRAM to fastq (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to CRAM, to FastQ (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    input:
        tuple val(sampleName), file(cram)

    output:
        tuple val(sampleName), path("${sampleName}_1.fastq.gz"), path("${sampleName}_2.fastq.gz")

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${sampleName}_1.fastq.gz -2 ${sampleName}_2.fastq.gz tmp.bam
        rm tmp.bam
        """
}

