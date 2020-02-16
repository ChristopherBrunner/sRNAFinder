#!/usr/bin/env nextflow

/*
* This sRNA window pipeline divides a whole bacterial genome into windows, that have potential to be small non coding sRNA.
* These windows are then analyzed by the sRNA-Characterization-Pipeline (sRNACharP) and ranked by the sRNA-Ranking script (sRNARanking).
* The output of these pipelines is then plotted for easier analyzation and interpretation.
*
* Input:
*   - Genome FASTA file, a GFF file, an organism name, the work directory for sRNACharP, a variable window size and step
*   - Static inputs are the RandomForest Classifier R script and rds file for sRNARanking, the sRNA Characterization Pipeline Nextflow script and the config file used to run it
*   - Optional inputs include the reference file, used for graphical comparision of the results, and the chromosome length file, splitting the genome into chromosome regions to be analyzed
* 
* Output:
*   - Two HTML files containing the plotted sRNA windows and their respective scores for being an sRNA or not being an sRNA
*   - The file used for plotting, containing the scores, strand, and genome position for each putative sRNA
*   - Genome window BED file containing the filtered sRNA window locations used to rerun the pipeline. This file can be used as the optional input file chromosome length,
*     to analyze genomic regions, that have a score above a certain threshold, with smaller windows in order to find shorter sRNAs.
*
* Author: Christopher Brunner
*/

// input files
params.org = ""                                         // Organism name
params.dir = ""                                         // Work directory for sRNACharP currently /home/christopher/Desktop/project-data/

if (params.org == "") {
  exit 1, "REQUIRED: Input an organism name as specified!"
}

if (params.dir == "") {
  exit 1, "REQUIRED: Input a work directory for sRNACharp containing the files as specified!"
}


params.fFile = "${params.org}_genome.fasta"             // Genome FASTA file
params.gffFile = "${params.org}_genome.gff"             // GFF file
params.RFClassifierR = "RF_classifier4sRNA.R"           // RandomForest R script
params.RFClassifierRDS = "RF_classifier4sRNA.rds"       // RandomForest rds file
params.sRNACharP = "sRNACharP.nf"                       // sRNACharP Nextflow code
params.sRNACharPConfig = "nextflow_sRNACharP.config"    // sRNACharP config file
params.refFile = ""                                     // Wet lab reference file
params.genomeLength = null                              // Optional genome length file, needed for a pipeline rerun or sRNA analyzation in one or multiple genomic regions
params.windowSize = 1000                                // Window size variable; standard set to 1000 nucleotides
params.step = 0.5                                       // Step variable is a fraction of the window size, which determines the distance between the windows; standard set to 0.5

fastaFile = file(params.fFile)
gffFile = file(params.gffFile)
RF_classifier4sRNA_R = file(params.RFClassifierR)
RF_classifier4sRNA_rds = file(params.RFClassifierRDS)
sRNACharP_nf = file(params.sRNACharP)
nextflow_sRNACharP_config = file(params.sRNACharPConfig)
windowSize = params.windowSize
step = params.step

chromLen = file("${params.org}ChromosomeLength.txt")

// Create a file with genomic start and end positions for each chromosome and split into two channels
Channel
    .from(fastaFile)
    .splitFasta( record: [id: true, seqString: true ])
    .collectFile(name: chromLen) { record -> record.id + "\t" + "0" + "\t" + record.seqString.length() + "\n"}
    .set{lengthGenome}

lengthGenome.into {
    lengthGenome_sRNAbed
    lengthGenomePlot
    lengthGenomeFilter
}

// In case a genome length file is given as input, overwrite existing file in the channel to cover only genomic regions specified in optional file
if (params.genomeLength != null) {
    genLen = file(params.genomeLength)

    Channel
    .from(genLen)
    .set{lengthGenome_sRNAbed}
}

/*
* Python script that creates a BED file containing the putative sRNA windows
* These entries are formatted as following: geneID, start, end, sRNA-name, score, strand
*/ 
process createsRNAbedFile {
    input:
        file "${params.org}GenomeLength.txt" from lengthGenome_sRNAbed
        val windowSize
        val step
    output:
        file "genomeSRNAs.bed" into sRNAWindowsBED

    """
#!/usr/bin/env python3

nameCount = 1

with open("${params.org}GenomeLength.txt", "r") as chromLen:
    with open("genomeSRNAs.bed", "w") as genomeSRNAs:
        # Function for adding one window to the sRNA BED file
        def printSRNA(geneID, start, end, strand):
            global nameCount
            genomeSRNAs.write(geneID + "\\t" + str(start) + "\\t" + str(end) + "\\tsRNA" + str(f"{nameCount:07}") + "\\t.\\t" + strand + "\\n")
            nameCount += 1

        # Iterate through all lines, each line contains one chromosome region that is supposed to be analyzed
        for gene in chromLen:
            geneID = gene.split()[0]
            start = int(gene.split()[1])
            end = int(gene.split()[2])
            maxStartPos = end - ${windowSize}

            # For each genomic region append an sRNA window to the BED file
            # Therefore iterate through the region with window size times the step to get a fraction of window size, creating a sliding window across the genome
            for index in range(start, end, int(${windowSize}*${step})):
                # Once the index exceeds the maximum start position, add one more window until the end of the region and then break out of the loop for this chromosome location
                if index > maxStartPos:
                    printSRNA(geneID, index, end, "+")  # Pos. strand
                    printSRNA(geneID, index, end, "-")  # Neg. strand
                    break
                else:
                    printSRNA(geneID, index, index+${windowSize}, "+") # Pos. strand
                    printSRNA(geneID, index, index+${windowSize}, "-") # Neg. strand
        """
}

// Split channel into two, one for running sRNACharP, the other for joining the sRNARanking result with their respective sRNA windows
sRNAWindowsBED.into {
    sRNAWindowsBEDsRNACharP
    sRNAWindowsBEDJoin
}

// Collect sRNA window BED file in working directory neccesary for running sRNACharP
sRNAWindowsBEDsRNACharP
.collectFile(name: file("${params.org}_genomesRNAWindows.bed"))

//Create the protein coding genome annotation file from the GFF file
process createProteinBEDFile {
    input:
         file gffFile
     output:
         file "genomeAnnotation.bed" into proteinCodingBED
     """
     awk 'OFS="\t" {if (\$3=="gene") {print \$1,\$4-1,\$5,\$9,".",\$7,\$3}}' ${gffFile}| perl -pe 's/\tID=.*Name=/\t/' | perl -pe 's/;.*\\w\t/\t/' > genomeAnnotation.bed
     """
 }

Collect genome annotation file in working directory neccesary for running sRNACharP
 proteinCodingBED
 .collectFile(name: file("${params.org}_genomeAnnotation_proteincoding.bed"))

// Run the sRNA characterization pipeline
process runsRNACharP {
    input:
        file sRNACharP_nf
        file nextflow_sRNACharP_config
    output:
        file "${params.org}_FeatureTable.tsv" into sRNACharPresult
    """
    nextflow run ${sRNACharP_nf} -c ${nextflow_sRNACharP_config} --org="${params.org}" --dir="${params.dir}" --genome="${params.org}_genome.fasta"  --sRNAs="${params.org}_genomesRNAWindows.bed" --genomeAnnotation="${params.org}_genomeAnnotation_proteincoding.bed"
    """
}

// Collect Feature Table in working directory neccesary for running sRNARanking
sRNACharPresult
.collectFile(name: file("${params.org}_FeatureTable.tsv"))

// Run the sRNA ranking R script
process runsRNARanking {
    input:
        file "${params.org}_FeatureTable.tsv" from sRNACharPresult
        file RF_classifier4sRNA_R
        file RF_classifier4sRNA_rds
    output:
        file "PipelineOutFile" into sRNARankingResult
    """
    Rscript --vanilla ${RF_classifier4sRNA_R} -i ${params.org}_FeatureTable.tsv -o PipelineOutFile
    """
}

/*
* Join the sRNA Ranking scores with their respective sRNA window's location, needed for the plot
* Add column headers for the plot
*/
process joinScoreAndCoordinates {
    input:
        file "PipelineOutFile" from sRNARankingResult
        file "genomeSRNAs.bed" from sRNAWindowsBEDJoin
    output:
        file "PipelineOutFileJOINED" into sRNARankingJoined
    """
    sort -k1 PipelineOutFile | perl -p -e 's/"//g' > PipelineOutFileSORTED
    join -1 1 -2 4 PipelineOutFileSORTED genomeSRNAs.bed | sed 's/ /\t/g' | sort -k1 > PipelineOutFileJOINED
    sed -i '1i\
    name\tNosRNA\tsRNA\tgeneID\tstart\tend\tscore\tstrand
    ' PipelineOutFileJOINED
    """
}

// Output channel split, one for filtering for the rerun of the pipeline, one for plotting and one as an output of the scores and regions in table format
sRNARankingJoined.into {
    sRNARankingJoinedFilter
    sRNARankingJoinedPlot
    sRNARankingJoinedOutput
}

sRNARankingJoinedOutput
.collectFile(name: file("${params.org}_Genome_sRNAScores.txt"))

/*
* Plot the scores for each sRNA window in a bar chart using csv, plotly and plotly express
* Each window is represented as one bar with a the height between 0 and 1 determening the score
* Two figures(files), one for the scores representing the windows being an actual sRNA (sRNA Values), the other representing the windows not being an actual sRNA (NosRNA Values)
* Each figure split into two subplots, respectively for each strand
*/
process plotGenomesRNAs {
    input:
        file "GenomesRNAScores" from sRNARankingJoinedPlot
        file "chromLen" from lengthGenomePlot
        val windowSize
        val step
    output:
        file "figuresRNA.html" into GenomePlotsRNA
        file "figureNosRNA.html" into GenomePlotNosRNA
    """
#!/usr/bin/env python3
import csv
import plotly.express as px
import plotly.graph_objects as go

sRNA_list = []

# Read lines of the joined ranking score and window location with csv
with open("GenomesRNAScores", "r") as sRNA_data:
    sRNA_reader = csv.DictReader(sRNA_data, delimiter="\t")
    
    # Append lines to sRNA_list as a dictionary
    for sRNA in sRNA_reader:
        sRNA_list.append(dict(sRNA))

# Change chromosome regions to genome coordinates
with open("chromLen", "r") as chromLen:   
    chromLenPrev = 0
    
    for chrom in chromLen:
        for sRNA in sRNA_list:
            if (sRNA['geneID'] == chrom.split()[0]):
                sRNA['start'] = int(sRNA['start']) + chromLenPrev
                sRNA['end'] = int(sRNA['end']) + chromLenPrev
        
        chromLenPrev = chromLenPrev + int(chrom.split()[2])

# Create a plotly express bar chart, one for sRNA Values and the other for NosRNA Values:
# sRNA Values
figuresRNA = px.bar(sRNA_list,                                      # List containing the dictionary for each sRNA window
                facet_row="strand",                                 # Splits each figure into two subplots, depending on the strand
                x="start",                                          # x-axis set to start poisiton
                y="sRNA",                                           # y-axis set to sRNA (value for being an sRNA)
                hover_data=["end", "name", "start", "strand"],      # Additional information added as hover data (end position, start position, name of the sRNA and the strand)
                labels={"sRNA":"Probability for an sRNA",           # Rename the x- and y-axis labels aswell as the end position
                 "start":"Start position", "end":"End position"}
)

# NosRNA Values
figureNosRNA = px.bar(sRNA_list,
                facet_row="strand",
                x="start", 
                y="NosRNA",                                              # y-axis set to NosRNA (value for not being an sRNA)
                hover_data=["end", "name", "start", "strand"], 
                labels={"NosRNA":"Probability for not being an sRNA",
                 "start":"Start position", "end":"End position"}          
)

# Read the reference data like the scores table with csv
if "${params.refFile}" != "":
    with open("${params.dir}/${params.refFile}") as ref_data:
        ref_reader = csv.DictReader(ref_data, delimiter="\t")
        
        # Add each reference sRNA
        for ref_sRNA in ref_reader:

            # Function to add the reference sRNA as a subplot to one of each of the figure's subplots
            def addTrace(figure, subplot):
                # Reference sRNA added as a red line on the x-axis
                figure.add_trace(go.Scatter(
                                        x=[ref_sRNA['start'], ref_sRNA['end']],
                                        y=[0, 0],
                                        hoverinfo="x",
                                        mode="lines",
                                        line=go.scatter.Line(color="red", width=15),
                                        showlegend=False),
                                        subplot,
                                        1)
            # If reference sRNA is on the positive strand, add the current subplot trace to the second subplot (pos. strand) of each figure
            if ref_sRNA['strand'] == '+':
                addTrace(figuresRNA, 2)
                addTrace(figureNosRNA, 2)
            # Otherwise add them to the first subplot of each figure (neg. strand)
            else:
                addTrace(figuresRNA, 1)
                addTrace(figureNosRNA, 1)

# Rename each of the figures title
figuresRNA.update_layout(title_text="sRNA Values for ${params.org} (window size: " + str(${windowSize}) + " and step size: " + str(int(${windowSize}*${step})) + ")")
figureNosRNA.update_layout(title_text="NosRNA Values for ${params.org} (window size: " + str(${windowSize}) + " and step size: " + str(int(${windowSize}*${step})) + ")")

# Fix each y-axis to the range 0 to 1
figuresRNA.update_yaxes(range=[0, 1])
figureNosRNA.update_yaxes(range=[0, 1])

figuresRNA.write_html("figuresRNA.html")
figureNosRNA.write_html("figureNosRNA.html")
    """
}

// Collect plot html files
GenomePlotsRNA
.collectFile(name: file("${params.org}_GenomePlotsRNA.html"))

GenomePlotNosRNA
.collectFile(name: file("${params.org}_GenomePlotNosRNA.html"))

/* 
* Trim the chromosome length file to only chromsome ID and length
* Filter out genome regions where the score for being a real sRNA is below the threshhold value of 0.405
* Extend each window up- and downstream by half of window size respectively (double window size in total) using bedtools slop
* Additionally filter out same region but different stranded windows
* Merge overlapping windows into one with bedtools merge
*/
process filterGenomeWindows {
    input:
        file "GenomesRNAScores" from sRNARankingJoinedFilter
        file "chromLen" from lengthGenomeFilter
        val windowSize
        val step
    output:
        file "FilteredGenomeWindows" into pipelineRerun
    """
    cut -f1,3 chromLen > onlyChromLen
    awk 'NR>1 && OFS="\t" {if (\$3 > 0.405) { print \$4, \$5, \$6 }}' GenomesRNAScores | uniq -f1 > FilteredGenome
    bedtools slop -i FilteredGenome -g onlyChromLen -b ${windowSize}/2 > ExtendedFilteredGenome
    bedtools merge -i ExtendedFilteredGenome > FilteredGenomeWindows
    """
}

// Collect Genome Windows file for the optional pipeline rerun
pipelineRerun
.collectFile(name: file("${params.org}_GenomeWindows.bed"))
