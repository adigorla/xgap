# XGAP

![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)
![release: v0.5.7](https://img.shields.io/badge/release-v0.5.7-green)
![coverage: 100%](https://img.shields.io/badge/coverage-100%25-brightgreen)
![docs: in-progress](https://img.shields.io/badge/docs-in--progress-yellow)

_**Documenation development in-progress**_

xGAP (e**x**tensible **G**enome **A**nalysis **P**ipeline) is an efficient, extensible, modular, open-source, fault-tolerant next generation sequencing(NGS) genome analysis pipeline framework. It implements massive parallelization of the GATK best practices pipeline  by splitting a genome into many smaller regions with efficient load-balancing to achieve high scalability. It is compatible for use on on-premises (SGE & SLURM) high performance clusters and AWS EC2 cloud servers (using [ParallelCluster](https://github.com/aws/aws-parallelcluster)). It is currently developed to efficiently perform germline small variants discovery accordance in with the [GATK best practices guidlines](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). xGAP fully automates each of processing steps in data analysis, from converting raw sequencing reads(FASTQs or BAMs) into variant call format (VCFs), and a optional VQSR step. This allows researchers with even limited computational backgrounds to take advantage of our pipeline. This pipeline allows rapidly generate germline small variant VCFs with minimal user intervantion and runtime. The end-to-end analysis runtime scales down as the number of core available (and number regions you can parallely process) increases. It is currently benchmarked to process 30x coverage whole-genome sequencing (WGS) data in approximately 90 minutes, with 256 cores. xGAP achieves average F1 scores of 99.37% for SNVs and 99.20% for indels across seven benchmark WGS datasets. To date xGAP has been validated for compatibility withe HG19, GRCh37 and GRCh38. The pipeline is easy-to-use and fault-tolerant because it can automatically re-initiate failed processes to minimize required user intervention without halting analyis. To use xGAP, one simply needs to update a runtime resource configuration file, provide paths to input files and execute xGAP using a simple command, after which the pipeline will fault-tolerantly progress without any further user intervention, sending a notification email after completion. The pipeline allows batch processing aswell. This allows you submit multiple WGS for analysis at the same time. The full software and source code is available to ensure anyone can use and modify the program freely to suit their research needs. 


## Installation

### Required Dependencies

* Python3
* pysam
* GATK 3.6+ or 4+
* Java 1.8.X
* BWA
* sambamba (only tested on v0.6.6 because newer versions have had issues on Intel Xeon processors)
* Conda (OPTIONAL)

### The Process

**1.** Clone or [Download](https://github.com/Adigorla/xgap/archive/master.zip) the latest version from this repo. Please ONLY clone the master branah as this is only version we can guarentee is stable. Clone/download the pipeline into your `$HOME` (or preferred directory).

```bash
git clone --single-branch --branch master git@github.com:Adigorla/xgap.git
```
**2.** Assign path to your local Python3 interpreter. 

```bash
~/xgap/interpAssign.sh </path/to/xgap/root/dir> </path/to/python3.X/interpreter>
```

**NOTE**: You only need to do this once after dowload to make sure xGAP knows where to look to run the python code. 

### Running xGAP consistently across multiple platforms

**Fill this in ..... **

## Usage

**1. Update configuration file**

This is the file where you set runtime resource and dependency configurations. This where you set the paths to external dependencies (such as GATK, BWA, etc.),  reference files (such as GRCh 37, dbsnp vcf, etc.) and intermediate/output file directories. This also the file where you configure scheduler runtime settings such as runtime and memory limits. Finally, users also get to select the number of regions to parallelize over and number of times a process can fail before analysis is halted for user intervention. Note: as a general rule of thunb, we think it's most efficient to set `n-regions` (the number of regions to parallelize over) to 2X the number of cores available for analysis. A template configuration file is already provided at `~/xgap/config/temp_config.yml`. You simplly needs to edit the file using you favourite text editor.

**2. Update input file**
This is the text file where you set the paths to input files (FASTQs or BAMs) for analysis. An example input file in provide at `~/xgap/bin/input_example.txt`. For bacth analysis follow the same format where each line should refer to single genome. Each line of the input file should refer to a single genome in the following format depending on the input format:

```
### for input pair-end FASTQ files

<Perfered Name for WGS dataset> </path/to/R1/fastq>,</path/to/R2/fastq>


### for input BAM/CRAM file

<Perfered Name for WGS dataset> </path/to/BAM> 


### Note that xGAP supports glob patterns for many files, as follows:

NA12878 /path/to/NA12878-illumina/NA12878_Rep_8_fastq/NA12878-Rep-8_S8_L00*_R1_001.fastq,/path/to/NA12878-illumina/NA12878_Rep_8_fastq/NA12878-Rep-8_S8_L00*_R2_001.fastq

```
Note: xGAP can, for now, only accept uncompressed FASTQ files with file extension format as `.fastq`.


**3. Run**

Executing xGAP is simple... 

```bash
~/xgap/bin/xgap <run_name> </path/to/input/file> </path/to/config/file>
```

`run_name` is simply name of batch and/or single file run, and can be anything. If you have activated email notifications, you should recive an alert when analyis is completed. The final VCFs should be found in `</path/to/output/dir/set/in/config>/xgap_output/<Perfered Name for WGS dataset>/vcf`.

## License and Disclaimer

xGAP is publicly released under the AGPL-3.0 license (full license text found [here](https://github.com/Adigorla/xgap/blob/master/LICENSE)). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
* [GATK](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)
* [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
* [sambamba](https://github.com/biod/sambamba/blob/master/LICENSE)
