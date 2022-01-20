#######################
# Generate Bootstraps #
#######################

#Generate a bootstrapped alignment from the real alignment.
#The original file is also copied to the 0th bootstrap file, so that it can easily undergo the same steps as the bootstrap files (1-N).
rule generate_bootstrap_alignments:
    input:
        "Data/{family}.mafft"
    output:
        temp("Data/{family}.mafft_bs:{bootstrap_n}")
    shell:
        """
        set +o pipefail;
        if [ "{wildcards.bootstrap_n}" == "0" ]; then
        	cp Data/{wildcards.family}.mafft Data/{wildcards.family}.mafft_bs:0
        else
        	Codes/generate_bootstrap_alignments.py {input} {wildcards.bootstrap_n}
        fi
        """

#########################
# Select Two-copy Genes #
#########################

#Take only the two copy-genes.
rule extract_two_copy_genes:
    input:
        "Data/{family}.mafft_bs:{bootstrap_n}"
    output:
        temp("Data/{family}:{bootstrap_n}.tcg")
    shell:
        """
        set +o pipefail;
        Codes/extract_two_copy_genes.py {input} {output}
        if [ "{wildcards.bootstrap_n}" == "0" ]; then
        	rm -f Data/{wildcards.family}.mafft_bs:0
        fi
        """

##############################
# Count Patterns in Quartets #
##############################

#Count phylogenetically informative positions (i.e. AABB, ABBA, etc.)
rule count_sequence_patterns:
    input:
        "Data/{family}:{bootstrap_n}.tcg"
    output:
        temp("Data/{family}:{bootstrap_n}.pts")
    shell:
        "Codes/count_sequence_patterns.py {input} {output}"

################
# Gene Network #
################

#Convert pattern counts to F (fate similarity measure) and create a network.
#Note that the networks can all be saved by removing the "temp(" and ")"; it is required in the next rule and later, to calculate the magnitude of recurrence.
#It can also be readily imported in Cytoscape.
rule make_gene_network:
    input:
        "Data/{family}:{bootstrap_n}.pts"
    output:
        temp("Data/{family}:{bootstrap_n}.fnet")
    shell:
        "Codes/make_fate_network.py {input} {output}"

#################
# Cluster Fates #
#################

#Cluster the network with Markov Clustering algorithm via an R-script (requires the packages: " ").
#Current parameters: expansion = 5, inflation = 10, edge_threshold = 0.0 (no threshold).
#Genes that do no have any connections fall out of the clustering so they need to be added as separate clusters by a few lines of bash code (see beneath "shell").
rule cluster_fates:
    input:
        fate_network="Data/{family}:{bootstrap_n}.fnet",
        two_copy_genes="Data/{family}:{bootstrap_n}.tcg"
    output:
        temp("Data/{family}:{bootstrap_n}.fcl")
    shell:
        """
        set +o pipefail;
        Codes/cluster_fates.R {input.fate_network} markov temp_{wildcards.bootstrap_n}.fcl 5 10 0.0
        cp temp_{wildcards.bootstrap_n}.fcl {output}
        ALL_TCG_GENES=`cat {input.two_copy_genes} | grep ">" | cut -c2-`
        for GENE in $ALL_TCG_GENES;
        do
            CHECK_PRESENCE=`grep ${{GENE:0:10}} temp_{wildcards.bootstrap_n}.fcl | wc -l`
            if [ "$CHECK_PRESENCE" -lt "1" ]; then
                echo -en $GENE"\n" >> {output}
            fi
        done
        rm -f temp_{wildcards.bootstrap_n}.fcl
        """

####################
# Merge Bootstraps #
####################

#Cluster bootstraps based on jaccard scores (REF), i.e. clustering of clusters.
#Here the number of bootstraps can be set! (in the brackets of "range()").
rule cluster_bootstraps:
    input:
        original="Data/{family}:0.fcl",
        bootstraps=expand("Data/{{family}}:{bootstrap_n}.fcl", bootstrap_n=[str(x+1) for x in range(100)])
    output:
        "Data/{family}.fates"
    shell:
        "Codes/cluster_bootstraps.py {output} {input.original} {input.bootstraps}"

###############
# Calculate P #
###############

#Calculate the prevalence of recurrent sequence evolution (i.e P, see paper).
rule measure_recurrence_prevalence:
    input:
        indies="Data/{family}.indies",
        fates="Data/{family}.fates"
    output:
        "Data/{family}.score"
    shell:
        "Codes/measure_recurrence_prevalence.py {input.indies} {input.fates} {output}"

#################
# Calculate Z_F #
#################

#First find the recurrent pairs in the largest fate of our family (i.e. the one that defines the family's P-score.
rule find_recurrent_quartets:
    input:
        indies="Data/{family}.indies",
        fates="Data/{family}.fates",
        two_copy_genes="Data/{family}:0.tcg",
        score="Data/{family}.score"
    output:
        non_rec_pairs1=temp("Data/{family}_nonrecpairs1.txt"),
        non_rec_pairs2=temp("Data/{family}_nonrecpairs2.txt"),
        recurrent_pairs=temp("Data/{family}_recurrentpairs.txt")
    shell:
        "Codes/find_recurrent_quartets.sh {input.indies} {input.fates} {input.two_copy_genes} {input.score} {output.non_rec_pairs1} {output.non_rec_pairs2} {output.recurrent_pairs}"

#Now calculate the magnitude of recurrent sequence evolution (i.e Z_F, see paper).
rule measure_recurrence_magnitude:
    input:
        score="Data/{family}.score",
        fates="Data/{family}.fates",
        patterns="Data/{family}:0.pts",
        non_rec_pairs1="Data/{family}_nonrecpairs1.txt",
        non_rec_pairs2="Data/{family}_nonrecpairs2.txt",
        recurrent_pairs="Data/{family}_recurrentpairs.txt"
    output:
        "Data/{family}.zscore"
    shell:
        "Codes/measure_recurrence_magnitude.sh {wildcards.family} {input.score} {input.fates} {input.patterns} {input.non_rec_pairs1} {input.non_rec_pairs2} {input.recurrent_pairs} {output}"

###############################
# Identify recurrent patterns #
###############################

#Here you will generate a reordered alignment file and relating annotation files for visualization in JalView.
rule identify_recurrent_patterns:
	input:
		alignment="Data/{family}:0.tcg",
		indies="Data/{family}.indies",
		fates="Data/{family}.fates"
	output:
		position_colors="Output/{family}_poscol.jalview",
		groups="Output/{family}_groups.jalview",
		reordered_alignment="Output/{family}.aln",
		seq_logo="Output/{family}_seqlogo.png"
	shell:
		"Codes/identify_recurrent_patterns.py {input.alignment} {input.indies} {input.fates} {output.position_colors} {output.groups} {output.reordered_alignment} {output.seq_logo}"

######################
# Collect all output #
######################

rule bottom_rule:
    input:
        one="Data/{family}.score",
        two="Data/{family}.zscore",
        three="Output/{family}_poscol.jalview"
    output:
        "Output/{family}.out"
    shell:
        """
        set +o pipefail;
        cat {input.one} {input.two} > {output}
        rm -f {input.one} {input.two}
        head -1 {input.three} >> check.out
        rm -f check.out
        """
