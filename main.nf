#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                            P I P E L I N E
========================================================================================
 Cellular Genetics bulk-RNA-Seq analysis pipeline, Wellcome Sanger Institute
 #### Homepage / Documentation
 https://github.com/cellgeni/
 #### Authors
 Vladimir Kiselev @wikiselev <vk6@sanger.ac.uk>
 Stijn van Dongen <svd@sanger.ac.uk>
 Original development by SciLifeLabs
----------------------------------------------------------------------------------------
*/

// def helpMessage() {
//     log.info"""
//     =========================================
//      Bulk-RNA-Seq pipeline v${version}
//     =========================================
//     Usage:

//     The typical command for running the pipeline is as follows:

//     nextflow run cellgeni/RNAseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile farm3

//     Mandatory arguments:
//       --genome                      Name of iGenomes reference
//       -profile                      Hardware config to use. farm3 / farm4 / openstack / docker / aws

//     Strandedness:
//       --forward_stranded            The library is forward stranded
//       --reverse_stranded            The library is reverse stranded
//       --unstranded                  The default behaviour

//     References                      If not specified in the configuration file or you wish to overwrite any of the references.
//       --star_index                  Path to STAR index
//       --star_overhang               sjdbOverhang parameter for building a STAR index (has to be (read_length - 1))
//       --fasta                       Path to Fasta reference
//       --gtf                         Path to GTF file
//       --bed12                       Path to bed12 file
//       --saveAlignedIntermediates    Save the BAM files from the Aligment step  - not done by default

//     Other options:
//       --outdir                      The output directory where the results will be saved
//       --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
//       -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
//     """.stripIndent()
// }


/*
 * utils is shared between projects. Include it in the PATH so scripts are found.
 */
// env.PATH = "$baseDir/utils:$PATH"


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.samplefile = false
params.studyid = -1
params.fastqdir = false
params.outdir = './results'
params.fcextra = ""                          // feature counts extra parameters; currently for testing
params.singleend = false

// Configurable variables
params.scratch = false
params.runtag  = "NF"                        // use runtag as primary tag identifying the run; e.g. studyid
params.name = false
params.project = false
params.genome = 'GRCh38'
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
params.salmon_trans_gene = params.genome ? params.genomes[ params.genome ].salmon_trans_gene ?: false : false
params.star_overhang = '74'
params.dna = params.genome ? params.genomes[ params.genome ].dna ?: false : false
params.cdna = params.genome ? params.genomes[ params.genome ].cdna ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.splicesites = false
params.download_hisat2index = false
params.download_fasta = false
params.download_gtf = false
params.hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.saveAlignedIntermediates = false
params.biotypes_header= "$baseDir/assets/biotypes_header.txt"


/*
 * Create a channel for input sample ids
 */
sample_list = Channel.fromPath(params.samplefile)

if (params.studyid > 0) {
    ch_fastqs_dir = Channel.empty()
    process irods {
        tag "${samplename}"

        input: 
            val samplename from sample_list.flatMap{ it.readLines() }
        output: 
            set val(samplename), file('*.cram') optional true into ch_cram_files
        script:
        """
        bash -euo pipefail irods.sh -t ${params.studyid} -s ${samplename}
        """
    }
} else if (params.fastqdir) {
    ch_cram_files = Channel.empty()
    if (params.singleend) {
      process get_fastq_files_single {
          tag "$samplename"

          input:
              val samplename from sample_list.flatMap{ it.readLines() }
          output:
              set val(samplename), file("${samplename}.fastq.gz") optional true into ch_fastqs_dir
          script:
          """
          name=${params.fastqdir}/${samplename}.fastq.gz
          if [[ ! -e \$name ]]; then
            echo "Count file \$name not found"
            false
          else
            ln -s \$name .
          fi
          """
      }
    }
    else {
      process get_fastq_files {
          tag "${samplename}"

          input:
              val samplename from sample_list.flatMap{ it.readLines() }
          output:
              set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_dir
          script:
          """
          list=( \$(ls ${params.fastqdir}/${samplename}_{1,2}.fastq.gz) )
          if [[ 2 == \${#list[@]} ]]; then
            ln -s \${list[0]} .
            ln -s \${list[1]} .
          else
            echo "Count mismatch sample ${samplename} found (\${list[@]})"
            false
          fi
          """
      }
    }
} else {
  exit 1, "Need --fastqdir <dirname> or --studyid <ID> option"
}


process crams_to_fastq {
    tag "${samplename}"

    if (params.scratch) {
       scratch true
    }

    input: 
        set val(samplename), file(crams) from ch_cram_files
    output: 
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_cram
    script:

        // 0.7 factor below: see https://github.com/samtools/samtools/issues/494
        // This is not confirmed entirely just yet.
        // def avail_mem = task.memory == null ? '' : "${ sprintf "%.0f", 0.7 * ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    """
    samtools merge -@ ${task.cpus} -f ${samplename}.cram ${crams}

    # check that the size of the cram file is >0.5Mb
    minimumsize=500000
    actualsize=\$(wc -c <"${samplename}.cram")

    f1=${samplename}_1.fastq.gz
    f2=${samplename}_2.fastq.gz

    if [ \$actualsize -ge \$minimumsize ]; then
                              # -O {stdout} -u {no compression}
                              # -N {always append /1 and /2 to the read name}
                              # -F 0x900 (bit 1, 8, filter secondary and supplementary reads)
      samtools collate    \\
          -O -u           \\
          -@ ${task.cpus} \\
          ${samplename}.cram pfx-${samplename} | \\
      samtools fastq      \\
          -N              \\
          -F 0x900        \\
          -@ ${task.cpus} \\
          -1 \$f1 -2 \$f2 \\
          -
    fi
    """
}


ch_fastqs_cram
  .mix(ch_fastqs_dir)
  .set{ ch_reads }

process assemble {
        tag "$samplename"
        publishDir "${params.outdir}", mode: 'copy',
            saveAs: { filename ->
                if (filename ==~ /.*\.ReadsPerGene\.out\.tab/) "STARcounts/$filename"
                else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
                else params.saveAlignedIntermediates ? "STARbams/filename" : null
            }

        input:
        set val(samplename), file(reads) from ch_reads
        file index from star_index.collect()
        file gtf from gtf_star.collect()

        output:
        set val(samplename), file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.SJ.out.tab"
        file "*.Log.out"
        file "*.Log.final.out" into star_log
        file "*.ReadsPerGene.out.tab"

        script:
                  // TODO featurecounts resorts the BAM file; SortedByName is not a STAR option though.
                  // --outSAMunmapped Within: In case someone wants the BAM files.
        """

        f1=${samplename}_1.fastq.gz
        f2=${samplename}_2.fastq.gz
        tracer assemble $f1 $f2 random_cell_name

        // STAR --genomeDir $index \\
        //     --sjdbGTFfile $gtf \\
        //     --readFilesIn $reads --readFilesCommand zcat \\
        //     --runThreadN ${task.cpus} \\
        //     --twopassMode Basic \\
        //     --outWigType bedGraph \\
        //     --outSAMtype BAM SortedByCoordinate \\
        //     --outSAMunmapped Within \\
        //     --runDirPerm All_RWX \\
        //     --quantMode GeneCounts \\
        //     --outFileNamePrefix ${samplename}.
        """
    }