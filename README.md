# ONT Unaligned Methylation Calling Pipeline
***in progress***

Written by Sophie Hoffman

This pipeline includes code to take raw ONT reads (filtered to a specific region of interest) as input and output a bam file of methylation probabilities for each read
1. Start with raw ONT reads
2. Locate reads of interest from a fasta file and nanopolish index
3. Filter to reads of interest using `fast5_subset`
4. Call methylation using guppy
5. Generate bam files with methylation probabilities using `guppy2sam`
6. Evaluate chromosome-specific methylation patterns in R

## Conda and Loading Packages

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

If you have trouble installing fast5mod, try this: 
```
pip install cython 
python -m pip install fast5mod --no-binary :all:
```

## 2. Locate reads of interest
In the use case this pipeline was designed for, we only wanted to call methylation on certain regions of the genome. If you want unaligned methylation calls for all reads in all fast5 files, skip this step. 
We identified the read IDs of interest by pulling them from a fasta file [bash script](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/pull_read_names.sh)

Then, using a [nanopolish](https://github.com/jts/nanopolish) index, which is created when you run nanopolish and tells you which fast5 file each read is located in, we select only the fast5 files we want to evaluate based on readnames of interest in this [R script](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/read_isolation.R)

## 3. Filter to reads of interest using `fast5_subset`
To extract the reads of interest from the fast5 files of interest, run our [bash loop](https://github.com/shoffm/ONT_unaligned_meth_calling/blob/master/fast5_sub_3.sh) (or better yet, parallize) to run [fast5_subset](https://github.com/nanoporetech/ont_fast5_api#fast5_subset)

This will generate [fast5 files](https://medium.com/@shiansu/a-look-at-the-nanopore-fast5-format-f711999e2ff6) containing the reads of interest. You can visually explore these files using [h5utils](https://github.com/NanoComp/h5utils).

## 4. Call methylation using guppy 

[Rerio](https://github.com/nanoporetech/rerio) allows you to use guppy to run methylation calling. 
- First, select a [model](https://github.com/nanoporetech/rerio#use-and-description-of-models) from the chart that suits your basecalling interests and sequencing devices used. If you have reads that were sequenced using both MinION/GridION and PromethION, select a model that is suitable for both, as not all models are. 
- [Install Rerio](https://github.com/nanoporetech/rerio#installation) and follow the instructions to download model(s) of interest. 



Visit this guide to [snakemake setup](https://github.com/Snitkin-Lab-Umich/Snakemake_setup) (including conda) to get started. 


