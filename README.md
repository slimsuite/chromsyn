# ChromSyn: Chromosome-level synteny plotting using orthologous regions

ChromSyn is designed to compile a set of BUSCO runs with the same version and lineage into chromosome synteny plots. This is achieved by establishing blocks of synteny based on co-linear regions that share an identifier. Whilst ChromSyn is designed with BUSCO in mind, it is therefore fairly simple to use alternative sources of synteny. Future releases will expand options for replacing BUSCO, so please get in touch if this would be useful to you.

## Version

The current version should be `v1.1.1`. (Check the chromsyn.R file to be sure!)

## Citation

If you use ChromSyn in your work, please cite the bioRxiv paper describing the Myrtle Rust genome: 

* Edwards RJ, Dong C, Park RF & Tobias PA Tobias. "A phased chromosome-level genome and full mitochondrial sequence for the dikaryotic myrtle rust pathogen, _Austropuccinia psidii_". bioRxiv 2022.04.22.489119 doi: [10.1101/2022.04.22.489119](https://doi.org/10.1101/2022.04.22.489119)

## Dependencies

You will need to have R with [tidyverse](https://www.tidyverse.org/) installed. The `writexl` package is also needed for Excel output.

## Quick start

The simplest way to use ChromSyn is using BUSCO predictions. Two inputs are required:

* `busco.fofn` = file of file names pointing to the BUSCO `full_table.tsv` tables
* `sequences.fofn` = file of file names pointing to sequence lengths

These are connected via the genome name for each assembly, and have contents in the simple form genome file with one line per genome:

```
GENOME1 FILENAME1
GENOME2 FILENAME2
...
GENOMEn FILENAMEn
```

One easy way to generate this input:

**Step 1.** Compile BUSCO results into a single directory and name them `$GENOME.busco5.tsv`.
	
**Step 2.** Name each input `$GENOME.fasta` and filter to just the chromosomes. (Can be done before or after running BUSCO but don't change sequence names if done afterwards.)
	
**Step 3.** Run [seqsuite](https://github.com/slimsuite/SLiMSuite) to build a sequences table per genome and copy the ``*.sequences.tdt` into the directory:

```
python $SLIMSUITE/tools/seqsuite.py -seqin $FASTA -seqmode db -summarise dna -basefile ${FASTA/.fasta/} i=-1
```

(Better still, run [Diploidocus](https://github.com/slimsuite/diploidocus) in telomere mode for telomeres to be plotted too - details and updated example output to follow.)

Version 0.6.0 introduced a new optional input of TIDK search output:

* `tidk.fofn` = file of file names for TIDK search results (`search/*_telomeric_repeat_windows.csv`)

Version 0.7.0 introduced a new optional input of assembly gaps (SeqSuite output) and features of interest:

* `gaps.fofn` = file of file names for assembly gaps (`*.gaps.tdt`)
* `ft.fofn` =  file of file names for features of interest (SeqName, Pos, [Strand,] [Col,] [Shape,])

If the features have one of the open shapes (21-25), they will have a black border and use the `Col` value for the fill.

**Step 4.** Make the FOFN files, e.g.:

```
for TSV in *.tsv; do
  echo ${TSV/.busco5.tsv/} $TSV >> busco.fofn
  echo ${TSV/.busco5.tsv/} ${TSV/.busco5.tsv/.sequences.tdt} >> sequences.fofn
done
```

**Step 5.** Edit `sequences.fofn` if required to set the vertical ordering of the species. (By default, the vertical ordering will match this file.)

**Step 6.** Run the chromsyn script:

```
Rscript $PATH/chromsyn.R basefile=$OUTPUTPREFIX | tee $OUTPREFIX.runlog
```

**Step 7.** Open up `$OUTPREFIX.pdf` to see the plot.


## Outputs

The main output is one or more [chromosome synteny plots](https://github.com/slimsuite/chromsyn/blob/main/zoomarsupials.pdf) and (if `writexl` is installed) an Excel file containing the established synteny. 

![Example plot of DNA Zoo marsupials. Each row is a different genome assembly. Each rectangle is an assembly scaffold. Reversed scaffolds have an R suffix. Default units are Mb, with tick marks every 10 Mb. Black and blue dots indicate Diploidocus and 3' TIDK telomere predictions. Other symbols mark annotated features, which are rRNA gene predictions in this example. Vertical placement indicates strand.](https://github.com/slimsuite/chromsyn/blob/main/zoomarsupials.png)

**NOTE:** ChromSyn will try to generate a single synteny plot. However, where there are many synteny blocks and rearrangements, system resources for R are sometimes exceeded. In this case, multiple plots will be produced that each contain a subset of the pairwise synteny comparisons. These can then be manually combined into a single figure.

**NOTE:** if provided, assembly gaps will be plotted as dark red plus signs.

## Options

A more detailed description of options and use cases will be added in time. The following options can be provided in the form `argument=value` to alter inputs and/or outputs:

```
# : sequences=FOFN = File of PREFIX FILE with sequence names and lengths (name & length, or SeqName & SeqLen fields) [sequences.fofn]
# : busco=FOFN = File of PREFIX FILE with full BUSCO table results. Used to identify orthologous regions. [busco.fofn]
# : tidk=FOFN = Optional file of PREFIX FILE with TIDK search results. [tidk.fofn]
# : regdata=TSV = File of Genome, HitGenome, Seqname, Start, End, Strand, Hit, HitStart, HitEnd
# : focus=X = If given will orient all chromosomes to this assembly
# : orient=X = Mode for sequence orientation (none/focus/auto)
# : seqsort=none/focus/auto/FILE = Optional ordering strategy for other assemblies [auto]
# : seqorder=LIST = Optional ordering of the chromsomes for the focal assembly
# : order=LIST = File containing the Prefixes to include in vertical order. If missing will use sequences=FOFN.
# : chromfill=X = Sequences table field to use for setting the colouring of chromosomes (e.g. Genome, Type or Col) [Genome]
# : basefile=FILE = Prefix for outputs [chromsyn]
# : plotdir=PATH = output path for graphics
# : minlen=INT = minimum length for a chromosome/scaffold to be included in synteny blocks/plots [0]
# : minregion=INT = minimum length for mapped regions to be included in plots [50000]
# : minbusco=INT = minimum number of BUSCO genes to be included in Syntenic block [1]
# : maxskip=0 = maximum number of BUSCO genes to skip and still be a syntenic block [0]
# : orphans=T/F = whether to include scaffolds that have no BUSCO genes [True]
# : tidkcutoff=INT = TIDK count cutoff for identifying a telomere [50]
# : align=X = alignment strategy for plotting chromosomes (left/right/centre/justify) [justify]
# : ygap=INT = vertical gap between chromosomes [4]
# : ypad=NUM = proportion of ygap to extend synteny blocks before linking [0.1]
# : scale=X = units in basepairs for setting the x-axis scale (bp/kb/Mb/Gb/pc) [Mb]
# : textshift=NUM = offset for printing chromosome names [0.3]
# : ticks=INT = distance between tickmarks [5e7]
# : pdfwidth=NUM = PDF width [20]
# : pdfheight=NUM = over-ride for standard calculated PDF height [0]
# : pdfscale=NUM = over-ride for PDF output scale [1]
# : ftsize=NUM setting to control the size of telomere and feature points.
# : namesize=NUM = scaling factor for the Genome names in PDF plots [1]
# : labelsize=NUM = scaling factor for the chromosome names in PDF plots [1]
# : labels=T/F = whether to print chromosome name labels [TRUE]
# : opacity=NUM = Opacity of synteny regions (0-1) [0.3]
# : debug=T/F = whether to switch on additional debugging outputs [FALSE]
# : dev=T/F = whether to switch on dev mode during code updates [FALSE]
```

## Running within R(Studio)

To use the script within other R code, the commandline arguments can be over-ridden with a vector of commandline arguments named `override`. For example, the code that generated the example plot:

```
override <- c("basefile=zoomarsupials","busco=zoomarsupials.busco.fofn","sequences=zoomarsupials.sequences.fofn","tidk=zoomarsupials.tidk.fofn","ft=zoomarsupials.ft.fofn","gaps=FALSE","orphans=F","minlen=1e7","focus=Wombat")
source("chromsyn.R")
```

## Algorithm overview

Chromosome synteny analysis is performed using single-copy ???Complete??? genes from BUSCO genome completeness runs. Synteny blocks between pairs of assembly are determined as collinear runs of matching BUSCO gene identifications. For each assembly pair, BUSCO genes rated as ???Complete??? in both genomes are ordered and oriented along each chromosome. Synteny blocks are then established as sets of BUSCO genes that are collinear and uninterrupted (i.e. sharing the same order and strand) in both genomes. A synteny block starts at the start of a BUSCO gene and will keep extending all the time the same BUSCO gene is found next in the same direction in both assemblies. It will end at the end of the last BUSCO gene in the block. (Local rearrangements and breakdowns of synteny between BUSCO genes will not be identified and may be falsely marked as syntenic within a synteny block.) 

Settings to control this behaviour include:

*	The min no. of genes needed to define a block. (1 by default)
*	The min length for a synteny block (50kb by default).
*	The max number of BUSCO genes that can be skipped to maintain synteny. (0 by default)


Chromosome synteny is then visualised by arranging the chromosomes for each assembly in rows and plotting the synteny blocks between adjacent assemblies. In each case, the chromosome order and orientation is set for one ???focus??? assembly, and the remaining assemblies arranged in order to maximise the clarity of the synteny plot, propagating out from the focal genome. For each chromosome, the ???best hit??? in each other assembly is established as that with the maximum total length of shared synteny blocks and an anchor point established as the mean position along the best hit chromosome of those synteny blocks. Where the majority of synteny is on the opposite strand, the chromosome is reversed and given an ???R??? suffix in the plot. Chromosomes are then ordered according to their best hits and, within best hits, the anchor points. Synteny blocks sharing the same orientation are plotted blue, whilst inversions are plotted in red.

 
