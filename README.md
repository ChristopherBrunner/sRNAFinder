# sRNAFinder

A pipeline for genome wide small RNA Characterization (sRNACharP) and Ranking (sRNARanking), written in the [Nextflow DSL](http://nextflow.io).

The pipeline calculates the scores for putative sRNA windows being a bona fide sRNA or not an sRNA at all, either on a genome wide scale or limited to one or multiple genomic regions. The output is presented as two plots, displaying the sRNA windows and their respective scores as bar graphs.

## 1. Getting Started

Since this pipeline is based on the two bioionformatics tools [sRNACharP](https://github.com/BioinformaticsLabAtMUN/sRNACharP) and [sRNARanking](https://github.com/BioinformaticsLabAtMUN/sRNARanking), install both tools first by following the instructions specified.
Make sure you download and install all of the files to the same directory in which you want to run this pipeline.
Keep in mind that remote runs on a high performance computing cluster are usually neccesary for genome wide analysis. Therefore Singularity will be used instead of Docker. (See 1.3)

### 1.1 Prerequisites

What you need to install in addition to the software needed for the bioinformatics tools mentioned above.

* [Python](https://www.python.org) 3.6 or higher
* [Singularity](https://sylabs.io/docs/) as a module on your server

### 1.2 Bprom issue

The bioinformatics tool bprom does not work when used from the container and needs to be used natively. Therefore, download the tool from [here](http://www.softberry.com/berry.phtml?topic=fdp.htm&no_menu=on) and add the direct path to the bprom file before the "bprom" command in the sRNACharP.nf process "getPromoterSites".

Additionally, the config file nextflow_sRNACharP.config needs to be adjusted with the path to the bprom/lin/data file.

Lastly, bprom will be replaced with a different tool that is still in development right now. This is, amongst other things, due to bprom producing an extremely high amount of files when used for multiple putative sRNAs, which may cause problems.

### 1.3 Docker and Singularity

The pipeline sRNACharP uses a Docker container, which can be used with Docker locally. Small genomic region runs of the sRNAFinder pipeline can be done locally with Docker. If you want to use Docker, move the config files from the Docker folder into the main working directory.  Make sure to delete the nextflow.config file provided by sRNACharP beforehand.

However, as mentioned above, most runs are very time consuming and and use a lot of memory. This is why it is recommended to use a high performance computing server, especially for whole genome runs. In this case Docker is replaced with Singularity.
Move the config files from the Singularity folder into the main working directory. Also make sure the the existing nextflow.config file is replaced with the new one.

Pulling the docker image file with Singularity on remote:

```
singularity pull docker://penacastillolab/srnacharp@sha256:c2a07f176d7cfe8cea3530bd76da05b30b182cdfe4d4b878f7d90e81f2d6a5f3
```
When using singularity, update both config files with the paths to the docker image where necessary.

## 2. Running the pipeline

### 2.1 Pipeline input
```
Options:
--org			Organism name (REQUIRED)
--dir			Directory containing the input files (REQUIRED)
--fFile			Genome FASTA file (REQUIRED)
--gffFile		GFF File (REQUIRED)
--RFClassifierR		RandomForest R script (REQUIRED)
--RFClassifierRDS	RandomForest rds file (REQUIRED)
--sRNACharP		sRNACharP Nextflow file (REQUIRED)
--sRNACharPConfig	sRNACharP Nextflow config file (REQUIRED)
--refFile		Wet lab reference file (OPTIONAL)
--chromosomeLength	Chromosome length file; can be used for pipeline re run or sRNA analyzation in one or multiple genomic regions (OPTIONAL: If not provided chromosome length file for whole genome will be calculated)
--windowSize		Window size variable; standard at 1000 nucletiodes (nts) (OPTIONAL)
--step			Step variable is a fraction of the window size, which determines the distance between the windows; standard set to 0.5 (OPTIONAL)
```
Make sure each input file is named correctly (```ORGANISM``` = organism name):

* fFile:             ```ORGANISM_genome.fasta```
* gffFile:           ```ORGANISM_genome.gff```
* RFClassifierR:	```RF_classifier4sRNA.R```
* RFClassifierRDS:	```RF_classifier4sRNA.rds```
* sRNACharP:		```sRNACharP.nf```
* sRNACharPConfig:	```nextflow_sRNACharP.config```

### 2.2 Pipeline run

The **default** pipeline run has a default window size of 1000 nts and a step of 0.5, meaning the pipeline will iterate through the genome with a 1000 nts sliding window that will step by a factor of 0.5 times the window size. In this example, this would mean a total step of 500 nts from the start of each window to the start of the next one and therefore yield an overlap between the sliding windows of 500 nts.

*(Ensure that the command is run in the correct working directory!)*

Here is the **default** pipeline run with the bacterium S. Enterica as an example:

```
nextflow run sRNAFinder.nf --org=”S_enterica” --dir=”/ADD_PATH/sRNAFinder/”
```
For **different values** of window size and step, use the following:
For example, use a step of 1 if no overlap between each sliding window of size 2000 is prefered.

```
nextflow run sRNAFinder.nf --org=”S_enterica” --dir=”/ADD_PATH/sRNAFinder/” --windowSize=”2000” --step=”1”
```
The pipeline is also capable of displaying **reference sRNAs** from the wet lab as a red bar on the x-axis of each plot.
To use this feature, include the following table headers as the first line of your reference BED file. Make sure that each column is named correctly and is tab-delimited. Keep in mind that each BED file can have varying numbers of columns, but start, end and strand have to be named exactly like in the example in order for this feature to work.

```
start    end    strand
```

Use this command to include the reference sRNAs in your plot, ”SLT2_sRNAs.bed” being the reference file name for S. Enterica

```
nextflow run sRNAFinder.nf --org=”S_enterica” --dir=”/ADD_PATH/sRNAFinder/” --refFile=”SLT2_sRNAs.bed”
```

After each run, the pipeline outputs a file named “ORGANISM_NAME_GenomeWindows.bed” with the genomic regions and gene ID of each putative sRNA that reached a threshold score for being a bona fide sRNA. These windows are extended in each direction by half their size (doubles each window size in total) and merged in case of an overlap between them.

This file can be used to **re run the pipeline** with a smaller window size, to detect smaller sRNAs that might be inside a larger sRNA window.

```
nextflow run sRNAFinder.nf --org=”S_enterica” --dir=”/ADD_PATH/sRNAFinder/” --windowSize=”250” --step=”1” --genomeLength=”S_enterica_GenomeWindows.bed”
```
The same input file can be used for **analyzing only one or multiple specific genomic regions**. Use the standard BED file format for this feature:

```
chromID    start	end
chromID	   start	end
chromID    start    	end
 ...        ...     	...
```
Finally, all options above can be used **together** or **seperately**, except for the **required** inputs –org and –dir

A **sample pipeline run** using all features would look something like the following:
```
nextflow run sRNAFinder.nf --org=”S_enterica” --dir=”/ADD_PATH/sRNAFinder/” --windowSize=”500” --step=”1” --genomeLength=”S_enterica_GenomeWindows.bed” --refFile=”SLT2_sRNAs.bed”
```

## 3. Output

* The pipeline output consists mainly of two HTML files containing the plotted putative sRNA windows and their respective scores, for being a bona fide sRNA or not being an sRNA
Plotly bar charts include hover data for each bar, containing detailed information of each sRNA window, a scrollable x axis and a selectable zoom
Note that in most cases a zoom in is needed to observe the bars, especially when plotting whole genomes
* Additionally, the file used for plotting, containing the scores, strand and genomic position of each sRNA window 

```ORGANISM_Genome_sRNAScores.txt```:
```
name	NosRNA	sRNA	geneID	start	end	score	strand
sRNA0000001	0.3775	0.6225	NC_003198.1	0	1000	.	+
sRNA0000002	0.685	0.315	NC_003198.1	0	1000	.	-
sRNA0000003	0.595	0.405	NC_003198.1	1000	2000	.	+
sRNA0000004	0.685	0.315	NC_003198.1	1000	2000	.	-
```
* A genome window bed file, containing the filtered sRNA window locations for the pipeline re run and is therefore used as the optional input chromosome length file

```ORGANISM_GenomeWindows.bed```:

```
NC_003198.1	0	2000
NC_003198.1	2000	4000
NC_003198.1	4000	6000
```

## Authors

* **Christopher Brunner**

## Acknowledgments

[sRNACharP](https://github.com/BioinformaticsLabAtMUN/sRNACharP) and [sRNARanking](https://github.com/BioinformaticsLabAtMUN/sRNARanking) created by Eppenhof EJ, Peña-Castillo L.
(2019) [Prioritizing bona fide bacterial small RNAs with machine learning classifiers. PeerJ 7:e6304](https://peerj.com/articles/6304/)
