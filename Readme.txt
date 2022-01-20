#README file to RecurrentEvolution workflow, by S.H.A. von der Dunk (2019); see the original paper.

The programme is currently a snakemake workflow comprised of several small scripts.
In the Data directory we included two input files as an example that allow you to run through the entire workflow by running "snakemake PTHR11822.out". Snakemake is a programme for building pipelines that works by making through specified rules the requested output-file (PTHR11822.out) if it can (i.e. given it has the input files PTHR11822.indies and PTHR11822.mafft). We refer to the Snakemake manuals for the details: https://snakemake.readthedocs.io/en/stable/



####################
### Installation ###
####################

In our case, the following steps had to be undertaken before the pipeline could be run on a private laptop with Ubuntu.
1) install miniconda using the website instructions.
2) install snakemake (with python 3.7) through conda.
3) install igraph through miniconda ("conda install -c r r-igraph")
4) install MCL from within R, installing it into the R-lib of miniconda: install.packages("MCL", "/home/sam/miniconda3/lib/R/library") [the latter path should obviously be substituted by the correct path in your own computer]



#################
### Execution ###
#################

See below for Installation.
Inside the main directory ("RecurrentEvolution") run the following command in the terminal (the input files for PTHR11822 are provided as an example):
> snakemake Output/PTHR11822.out

To have the whole workflow run snakemake on multiple families and with 6 cores:
> snakemake -j 6 Output/FAM1.out Output/FAM2.out Output/FAM3.out

Before running the actual pipeline it can be handy to check if snakemake can find the right input files; just add -np to the command you are going to run:
> snakemake -np -j 6 Output/FAM1.out Output/FAM2.out Output/FAM3.out

During the workflow, files always contain the family name in their name. Thus, one can in principle run the pipeline for multiple different families at the same time. However, we want to warn that snakemake tends to become REALLY inefficient when the number of jobs becomes very big. So, if you also like to do many bootstraps for each family, snakemake might become the biggest blockade (just run one family at a time, or call snakemake independently for each).
Normal execution can be performed by typing "snakemake PTHR11822.out" in the main directory.
By default, the number of bootstraps is set to 10. Change this in the Snakefile, under the "Merge Bootstraps" section, to do fewer or more bootstraps. Note that 1 bootstrap is the minimum, otherwise some rules do not get the right number of files for execution. Also note that the bootstraps bring in stochasticity. If different clusters emerge from the bootstrapped alignments, the scores (P and ZF) may also change.

To give an idea of the runtime of the whole workflow from scratch, we timed runs with the PTHR11822 example at a private laptop (excluding the last script):
1 bootstrap - 20s
10 bootstraps - 40s
100 bootstraps - 4m03s

Note that bootstraps generally barely alter the results (e.g. the P- and ZF-score). So when the workflow is used as a first-hand look into a family, one can probably just as well run with a few bootstraps only.



#############
### Input ###
#############

The input files serve as an example for usage of the pipeline.

The current workflow only implements the "clustering of fates" part, since duplications will generally be provided by the user through gene tree reconstruction.
The format of [family].indies seems a bit weird, but this stems from the fact that it used to be created by a workflow quite similar to the present pipeline.
Thus, to provide inferred duplications of a family, write all species that share a duplication on one line, separated by a tab, and preceeded by some number.
Also, the first line will be skipped by the scripts, so you have to write something (e.g. see the "Bootstrap = 100" in PTHR11822.indies).

The other input file is fairly straight-forward; it is the alignment of the family in FASTA format.



##############
### Output ###
##############

After running the example of PTHR11822 (see Execution), 5 files should appear in the Output directory:
	1) PTHR11822.aln: This contains the alignment of only the duplicates that underwent recurrent sequence evolution, grouped by duplications and ordered by fate (i.e. the first of two human paralogs (HSAP) has the same fate as the first of two Dictyostelium paralogs (DDIS)). Open this alignment in JalView.
	2) PTHR11822_groups.jalview:	This is one of the two annotation files for PTHR11822.aln in jalview. It draws boxes around duplicates that derive from the same duplication. Upload this file as annotation-file after opening PTHR11822.aln in JalView (note that not all JalView versions work well with annotation-files, in some cases restarting JalView can help).
	3) PTHR11822.out:	This shows the two largest predicted fate clusters (same as in Data/PTHR11822.fates) and the two family-level scores for recurrent sequence evolution (see paper); the number of independent duplications that fall completely in these two fates (P) and the average significance of the fate score (ZF) between duplicates in these fates that are not from the same duplication.
	4) PTHR11822_poscol.jalview:	This is the second of two annotation files for PTHR11822.aln in JalView. It colors each residue according to its consistency with the overall fate prediction (see paper).
	5) PTHR11822_seqlogo.png:	This is a sequence logo that summarizes the consistency of each position with the fate differentiation (i.e. integrating the different colors that PTHR11822_poscol.jalview assigns to residues for an entire column/position in the alignment).

Note that this output is also available for all the families that we analysed in the paper, so for these families it is not necessary to run the analysis again. Check out the Supplementary_Data.zip available from the same location as this present package.

