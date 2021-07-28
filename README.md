# ONT Unaligned Methylation Calling Pipeline
***in progress***

Written by Sophie Hoffman

This pipeline includes code to take raw ONT reads in the form of [fast5 files](https://medium.com/@shiansu/a-look-at-the-nanopore-fast5-format-f711999e2ff6) (filtered to a specific region of interest) as input and output a bam file of methylation probabilities for each read. We also include R scripts to process and visualize the resulting bam files into chromosome- and region-specific groups. 

1. Start with raw ONT reads, [load packages](#1-conda-and-loading-packages)
2. [Locate reads of interest from a fasta file and nanopolish index](#2-locate-reads-of-interest)
3. [Filter to reads of interest using `fast5_subset`](#3-filter-to-reads-of-interest-using-fast5_subset)
4. [Snakemake pipeline](#4-snakemake-pipeline)

    a. [Call methylation using guppy](#a-call-unaligned-methylation-using-guppy)
    
    b. [Generate bam files with methylation probabilities using `guppy2sam`](#b-generate-bam-files-with-methylation-probabilities-from-fast5-files-using-guppy2sam)
    
    c. [Merge bam files](#c-merging-bam-files)
    
    d. [Running our snakefile](#d-running-our-snakefile)
    
5. [Evaluate chromosome-specific methylation patterns in R](#5-evaluate-chromosome-specific-methylation-patterns-in-r)

## 1. Conda and Loading Packages

To create the conda environment:
```
conda env create meth_calling # you only have to do this once
conda activate meth_calling # you have to do this every time 
```
Next, install the following packages while in the conda environment
- [ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api#getting-started)
- [guppy](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html)
- [samtools](https://github.com/samtools/samtools)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [h5utils](https://github.com/NanoComp/h5utils) (for visualizing fast5 files) 
- [fast5mod](https://github.com/nanoporetech/fast5mod) 

If you have [trouble installing fast5mod](https://github.com/nanoporetech/fast5mod/issues/9#issuecomment-826920740), try this: 
```
pip install cython 
python -m pip install fast5mod --no-binary :all:
```

## 2. Locate reads of interest
In the use case this pipeline was designed for, we only wanted to call methylation on certain regions of the genome. If you want unaligned methylation calls for all reads in all fast5 files, skip this step. 
We identified the read IDs of interest by pulling them from a fasta file using this [bash script](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/pull_read_names.sh)

Then, using a [nanopolish](https://github.com/jts/nanopolish) index, which is created when you run nanopolish and tells you which fast5 file each read is located in, we select only the fast5 files we want to evaluate based on readnames of interest in this [R script](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/read_isolation.R)

We also generated a read reference file using this [bash script](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/pull_chm_spec_read_names.sh) so that we can later reference which region each read belongs to. 

## 3. Filter to reads of interest using `ont_fast5_api` [fast5_subset](https://github.com/nanoporetech/ont_fast5_api#fast5_subset)
To extract the reads of interest from the fast5 files of interest, run our [bash loop](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/fast5_sub_3.sh) (or better yet, parallize) to run [fast5_subset](https://github.com/nanoporetech/ont_fast5_api#fast5_subset)

This will generate [fast5 files](https://medium.com/@shiansu/a-look-at-the-nanopore-fast5-format-f711999e2ff6) containing the reads of interest. You can visually explore these files using [h5utils](https://github.com/NanoComp/h5utils).

## 4. Snakemake Pipeline

The following steps occur in our snakemake pipeline

### a. Call unaligned methylation using guppy

[Rerio](https://github.com/nanoporetech/rerio) allows you to use [guppy](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html) to run methylation calling on raw ONT reads (fast5 files). 
- First, select a [model](https://github.com/nanoporetech/rerio#use-and-description-of-models) from the chart that suits your basecalling interests and sequencing devices used. If you have reads that were sequenced using both MinION/GridION and PromethION, select a model that is suitable for both, as not all models are. See [navigating multiple sequencing platforms](#navigating-multiple-sequencing-platforms)
- [Install Rerio](https://github.com/nanoporetech/rerio#installation) and follow the instructions to download model(s) of interest. 

The guppy command we use in the snakefile (where `$` indicates a variable that you should change to suit your analysis):
```
guppy_basecaller --save_path $outpath_guppy --input_path $inpath_fast5_directory --compress_fastq --fast5_out --data_path $guppy_config_path --config $guppy_config; 
```

### b. Generate bam files with methylation probabilities from fast5 files using `guppy2sam`
`guppy2sam` allows you to convert a bam file from a fast5 file. 
The `guppy2sam` command we use in the snakefile (where `$` indicates a variable that you should change to suit your analysis): 
```
fast5mod guppy2sam $outpath_guppy/guppy_outfile.fast5 --workers 74 --recursive | samtools sort -@ 8 | samtools view -b -@ 8 > $outpath_bams;
```

### c. Merging bam files
In order to compile reads from multiple fast5 files, we merge the bam files using the following command in the snakefile (where you should change the `/path/to/bams/*` path): 
```
files=$(ls -d /path/to/bams/*)
samtools cat $files -o $outpath_merged_bams/all.bam
```

### d. Running our snakefile

Visit this guide to [snakemake setup](https://github.com/Snitkin-Lab-Umich/Snakemake_setup) to get started.

Edit our [snakefile_g](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/snakefile_g) file to include the locations of your input and desired output location in the designated locations. 
Run the following command as a dry run for snakemake, to ensure the paths you input work with the pipeline. 

```
snakemake -s snakefile_g -n
```

Run the snakemake pipeline as a job using the [snakemake.sh](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/snakemake.sh) submission script, or run it locally using the command in the file. 


## 5. Evaluate chromosome-specific methylation patterns in R
Using the newly generated bam file of all reads of interest, we can evlauate chromosome-specific methylation patterns of 45s rDNA units using this [R script](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/45_s_parse.R), where we find the probability of 5mC at CpG sites. We can visualize the results in this [R Markdown file](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/Chromosome-specific_methylation_analysis.Rmd)
*Note, this code subsets to regions within the reads from an external file, you can adjust it to define your own windows.*

## Bam file content explanation
The unaligned bam file generated will have: 
- A header
- Sequence
- Quality scores
- Methylation probabilities, tags for these sections will differ depending on which rerio model was used. I called 5mC, so the tag I wanted was MC:B:C. In this segment there are numbers between 0-255 which represent the probability that given a C was called, it was a methylated. 0 represents a 0% chance of 5mC, 255 represents a 100% chance of 5mC, and anything in between represents a (n/255)(100)%

## Navigating multiple sequencing platforms
Notes to handle reads from both MinION/GridION and PromethION: 

Ensure you choose a rerio model that exists for both devices, not all models have config files for both. 

In order to run guppy on the different set of reads, I was able to identify "promethion" in the filepath of the fast5 files that were run on promethion. 
I separated the filtered reads based on their original filepath into [separate directories](https://stackoverflow.com/questions/14371039/copy-all-files-in-a-directory-with-a-particular-string-in-the-filename-to-differ/14371380) and ran two different snakefiles, with a compatible rerio config file and directory of fast5 reads. *There is probably a better way to do this with a condition in the snakefile based on the fast5 filename, this was just the workaround I used.*

Acknowledgements: 
Initial inspiration came from this [aligned basecalling tutorial](https://labs.epi2me.io/notebooks/Modified_Base_Tutorial.html)

