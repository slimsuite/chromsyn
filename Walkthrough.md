# ChromSyn walkthrough

This walkthrough is designed to introduce the core functionality and primary use case for the ChromSyn synteny plotter. The [Tutorial](Tuorial.md) page will introduce some more use cases.

ChromSyn is designed to compile a set of BUSCO runs with the same version and lineage into chromosome synteny plots. This is achieved by establishing blocks of synteny based on co-linear regions that share an identifier. Whilst ChromSyn is designed with BUSCO in mind, it is therefore fairly simple to use alternative sources of synteny. 

## Installation and Setup

To run ChromSyn from scratch, you will need [Git](https://github.com/git-guides/install-git) and [R](https://cran.rstudio.com/) installed on your system with the the following R packages installed:

```
tidyverse
RColorBrewer
gtools
writexl
```

Then clone the [ChromSyn GitHub repository](https://github.com/slimsuite/chromsyn) and you are ready to run on existing data.

```
git clone https://github.com/slimsuite/chromsyn.git
Rscript chromsyn/chromsyn.R
```

If you do not have the required packages installed, the script should crash and tell you so. Otherwise, it will report you are missing required files:

```
[Wed Feb  8 15:15:15 2023] Sequence FOFN File: sequences.fofn
[Wed Feb  8 15:15:15 2023] Cannot find sequence FOFN file: sequences.fofn
```

If you already have data generated, or want to run on the example data provided, skip ahead to **Running ChromSyn on existing files**.

### Applications required for generating input files

In order to generate the full set of input files as described in this walkthrough, you will also need to have the following installed:

* [BUSCO](https://busco.ezlab.org/)
* [Diploidocus](https://github.com/slimsuite/diploidocus). (Clone the Git repo, no installation needed.)
* [TIDK](https://github.com/tolkit/telomeric-identifier)


## Generating input files

ChromSyn is designed to run from a set of chromosome-level genome assemblies in fasta format. Whilst the formatting and naming of these files is flexible, the recommended approach is to generate a fasta file of just those sequences of interest (e.g. the chromosome scaffolds) and rename them for clear labelling in the plots. For example, the data in the `example/gendata` folder was generated from three Ensembl marsupial reference genomes:

```
├── monodelphis_domestica
│   └── fasta
│       └── Monodelphis_domestica.ASM229v1.dna.toplevel.fa
├── ornithorhynchus_anatinus
│   └── fasta
│       └── Ornithorhynchus_anatinus.mOrnAna1.p.v1.dna.toplevel.fa
└── sarcophilus_harrisii
    └── fasta
        └── Sarcophilus_harrisii.mSarHar1.11.dna.toplevel.fa
```

These were then filtered to extract just the chromosome scaffolds (e.g. using name or length filtering) and renamed with simple chromosome numbering prefixes, e.g.:

```
sed -E 's/>(.+ chromosome )([0-9X])/>MONDOCHR\2 \1\2/' Monodelphis_domestica.ASM229v1.dna.toplevel.chrom.fa > ensMONDO.fasta 
```

In this case, this produced fasta files for three species:

* `ensDEVIL.fasta`: Sarcophilus harrisii = Tasmanian Devil. Chromosome names: `DEVILCHR1`, `DEVILCHR2`, `...`
* `ensMONDO.fasta`: Monodelphis domestica = Gray short-tailed opossum. Chromosome names: `MONDOCHR1`, `MONDOCHR2`, `...`
* `ensPLATY.fasta`: Ornithorhynchus anatinus = Duck-billed platypus. Chromosome names: `PLATYCHR1`, `PLATYCHR2`, `...`

### Generating sequence input

ChromSyn needs a file with the sequence names and lengths. This could be a simple `name,length` CSV file, but the recommended thing is to run [Diploidocus](https://github.com/slimsuite/diploidocus) and generate telomere predictions at the same time (where `$PATH` is the Dipoidocus code path). Diploidocus can also be used to create the gaps table:

```
for GENOME in *.fasta; do
  PREFIX=$(basename ${GENOME/.fasta/})
  python $PATH/diploidocus.py runmode=telomeres telonull seqin=$GENOME basefile=$PREFIX
  python $PATH/diploidocus.py runmode=summarise seqin=$GENOME gapstats dna i=-1 basefile=$PREFIX 
done
```

This should generate a `*.telomeres.tdt` file per genome, which will be used for the `sequences.fofn` input (see below), and `*.gaps.tdt` table for the `gaps.fofn` input.


### Running BUSCO

Next, we need the BUSCO synteny linkage. This will be run with a command similar to:

```
LINEAGE=/data/bio/busco/5/lineages/mammalia_odb10
PPN=40
for GENOME in *.fasta; do
  PREFIX=$(basename ${GENOME/.fasta/})
  busco -o run_$PREFIX -i $GENOME -l $LINEAGE --cpu $PPN -m genome
done
```

Obviously, your lineage path (and selection) will vary.

Once BUSCO has run, consolidate the `full_table.tsv` output from BUSCO:

```
for GENOME in *.fasta; do
  PREFIX=$(basename ${GENOME/.fasta/})
  cp -v run_$PREFIX/run_mammalia_odb10/full_table.tsv $PREFIX.busco5.tsv  
done
```

### Running TIDK

The final generic input for ChromSyn is the TIDK telomere repeat scores, in this case using a `AACCCT` telomere repeat:

```
for GENOME in ../*.chrom.fasta; do
  PREFIX=$(basename ${GENOME/.fasta/})
  tidk search -f $GENOME -o $PREFIX -s AACCCT
  cp -v search/${PREFIX}_telomeric_repeat_windows.csv $PREFIX.tidk.csv
done
```


## Running ChromSyn on existing files

Once the files have been generated (above), consolidate them in a directory as with the `example/gendata` ("genome data") directory provided:

```
example/gendata/
├── ensDEVIL.busco5.tsv
├── ensDEVIL.gaps.tdt
├── ensDEVIL.telomeres.tdt
├── ensDEVIL.tidk.csv
├── ensMONDO.busco5.tsv
├── ensMONDO.gaps.tdt
├── ensMONDO.telomeres.tdt
├── ensMONDO.tidk.csv
├── ensPLATY.busco5.tsv
├── ensPLATY.gaps.tdt
├── ensPLATY.telomeres.tdt
└── ensPLATY.tidk.csv
```

**NOTE:** This is not strictly necessary, as the `*.fofn` input files can point anywhere, but it makes their generation much easier.

Next, generate the `*.fofn` files for the four main input types. These all have the same simple structure:

```
GENOME1 FILENAME1
GENOME2 FILENAME2
...
GENOMEn FILENAMEn
```

If generated as above, with consistent file names that match the assembly names we want to use, the FOFN files can be made quite easily, e.g.:

```
cd gendata/
ls *.busco5.tsv | sed -E 's#([A-Za-z]+)\.#\1 gendata/\1.#' | tee ../busco.fofn
ls *.gaps.tdt | sed -E 's#([A-Za-z]+)\.#\1 gendata/\1.#' | tee ../gaps.fofn
ls *.telomeres.tdt | sed -E 's#([A-Za-z]+)\.#\1 gendata/\1.#' | tee ../sequences.fofn
ls *.tidk.csv | sed -E 's#([A-Za-z]+)\.#\1 gendata/\1.#' | tee ../tidk.fofn
```

In this case, the `*.fofn` files were created in the `example` directory and designed to be run from there, so the relative paths include the `gendata/` path:

```
ensDEVIL gendata/ensDEVIL.busco5.tsv
ensMONDO gendata/ensMONDO.busco5.tsv
ensPLATY gendata/ensPLATY.busco5.tsv
```

If running from another directory, replace `gendata` with the appropriate relative or absolute path. (See, for example, the `*.fofn` files in the `testrun/` directory.)

### Running with default settings

ChromSyn is now ready to run. To try out on the example data with defaults, enter the `testrun` data and run the R script:


```
Rscript ../chromsyn.R
```

If R was installed correctly and the github repo cloned OK, the following three files should have been created:

```
├── chromsyn.xlsx
├── chromsyn.pdf
└── chromsyn.png
```

This should produce an [example PDF plot](example/chromsyn.pdf) and the corresponding PNG:

![_Test plot. Each row is a different genome assembly. Each rectangle is an assembly scaffold. Reversed scaffolds have an R suffix. Default units are Mb, with tick marks every 10 Mb. Black and blue dots indicate Diploidocus and 3' TIDK telomere predictions. Vertical placement indicates strand._](example/chromsyn.png)


To capture the log output, redirect the stdout into a file:

```
Rscript ../chromsyn.R | tee chromsyn.log
```

This is not essential, but would be useful if things don't behave as expected.


