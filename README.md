### Introduction
---
A pipeline for identifying different derived neoantigens using RNA-Seq data and LC-MS data

### Usage
---
##### STEP 1: Data preprocessing and HLA typing acquisition.
`bash 00.Alignment/run.Alignment.hg38 sample fq1 fq2 outdir `

`bash 01.HLAtyping/run.OptiType-seq2HLA sample fq1 fq2 dir `

    -- sample: name of the sample
    -- fq1/fq2: the RNA sequencing fastq files
    -- outdir: path to output directory

##### STEP 2: Construction of cancer-specific database.
`bash 02.Cancer.specfic.proteome/kmerGeneration.sh sample fq1 fq2 outdir NB_THREADS`

    -- sample: name of the sample
    -- fq1/fq2: the RNA sequencing fastq files
    -- outdir: path to output directory
    -- NB_THREADS: number of processors to be used

`bash 02.Cancer.specfic.proteome/kmerFiltering.sh SAMPLE_NAME PATH_TO_JF LCOUNT PATH_TO_JF_FILTER PATH_OUT`

    -- SAMPLE_NAME: name of considered sample
    -- PATH_TO_JF: path to jf database of 33-nt-long cancer k-mers
    -- LCOUNT: minimal occurence for cancer k-mers to be considered expressed (≥), needs to be adjusted to get a maximum of 30M cancer-specific k-mers
    -- PATH_TO_JF_FILTER: list of paths to jf database of 33-nt-long normal k-mers. k-mer filtering will be sequential
    -- PATH_OUT: path to output directory

`bash 02.Cancer.specfic.proteome/kmerAssembly.sh PROJECT_FILE KMER_FILE PATH_OUT LCOUNT`

    -- PROJECT_FILE: path to project file 
    -- KMER_FILE: path to k-mer file 
    -- PATH_OUT: path to output directory
    -- LCOUNT: same as the previous step, as an identifier in the output file name

##### STEP 3: Construction of sample-specific database including different mutation types.

`bash 03.mutation.specific.proteome/indel-fasta.sh sample outdir bam`

`bash 03.mutation.specific.proteome/snp-fasta.sh sample outdir bam`

`bash 03.mutation.specific.proteome/AS-fasta.sh sample AS_type bam outdir`

`bash 03.mutation.specific.proteome/genefusion-fasta.sh sample outdir fq1 fq2`

`bash 03.mutation.specific.proteome/edit-fasta.sh sample outdir bam`

    -- sample: name of the sample
    -- outdir: path to output directory
    -- bam: The BAM file obtained from alignment in STEP 1
    -- AS_type: Different types of alternative splicing, such as A3SS
    -- fq1/fq2: the RNA sequencing fastq files

##### STEP 4: Identification of TSA candidates.
`bash 04.Identification/capaMHC.sh sample MS_filter_file hla outdir`

`bash 04.Identification/annotation.sh sample outdir AS_type LCOUNT`

    -- sample: name of the sample
    -- MS_filter_file: TXT files obtained through mass spectrometry filtering
    -- hla: HLA typing file in the following format:
        HLA-A24:02
        HLA-B68:01
        HLA-B35:03
        HLA-C04:01
    -- outdir: path to output directory
    -- AS_type: Different types of alternative splicing, such as A3SS
    -- LCOUNT: same as the STEP2, as an identifier in the file name
