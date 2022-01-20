#!/usr/bin/env Rscript
###############################################################################################
#Start-up: library, read data
#http://kateto.net/networks-r-igraph for very clear tutorial on igraph

# Run this script with: fate_clustering.R [network-file] [markov] [clusters-output-file] [par_expansion] [par_inflation] [par_threshold]
# For example: fate_clustering.R family_123.net family_123.clusters 5 10 0.0
# The first two parameters are only for MCl (so if option "markov" is given), but give them anyway to please R.
# For methods other than MCl, see below under Options.
# The script will output the clusters it found, one on each line, and the elements within each cluster are separated by a tab.
# The network input should be according to iGraph (I think):
# Each edge is on a separate line: node1 "tab" edge_value "tab" node2
###############################################################################################

suppressMessages(library(igraph))
suppressMessages(library(MCL))
rm(list=ls())
args=commandArgs(trailingOnly = T)
readnetwork=toString(args[1])
if(file.info(readnetwork)$size == 0){
	sink(args[3], append=T, split=F);
	cat("\n")
	sink();
}else{
	connections <- read.table(readnetwork, stringsAsFactors = FALSE)

	cluster_method=toString(args[2])

	par_expansion=strtoi(args[4])
	par_inflation=strtoi(args[5])
	par_threshold=as.double(args[6])

	#Options:
	### markov (mcl)
	### labelprop
	### highconn (hcs)
	### community-fastgreedy  ...See: https://www.sixhat.net/finding-communities-in-networks-with-r-and-igraph.html
	### community-spinglass  ...Only works on a connected graph.
	### community-edgebetweenness
	### community-walktrap
	### noclustering

	Edges <- data.frame(from=connections$V1, to=connections$V3, strength=connections$V2, stringsAsFactors = FALSE)
	Threshold <- par_threshold #Parameter for the clustering!

	###############################################################################
	if(cluster_method!="highconn"){
		Network <- graph_from_data_frame(Edges[Edges$strength>=Threshold,1:3], directed=FALSE)
		Net_Length <- length(V(Network)$name)
	}else{
		gV <- unique(c(connections$V1,connections$V3))
		ngV <- 1:length(gV)
		gE <- vector("list", length=length(gV))
		names(gE) <- gV
		for(p in 1:length(gV)){
			edge_labels <- c(connections$V3[connections$V1==gV[p]],connections$V1[connections$V3==gV[p]])
			edgesL <- rep(0,length(edge_labels))
			for(r in 1:length(edge_labels)){
				edgesL[r] <- ngV[gV==edge_labels[r]]
			}
			edges <- edgesL
			weights <- c(connections$V2[connections$V1==gV[p]],connections$V2[connections$V3==gV[p]])
			gE[[p]] <- list(edges=edges, weights=weights)
		}
		Network <- graph::graphNEL(nodes=gV, edgeL=gE)
		Net_Length <- length(gV)
	}
	#############################################################################

	if(Net_Length != length(unique(c(connections$V1, connections$V3)))){
		No_Lost_Vertices <- length(unique(c(connections$V1, connections$V3))) - Net_Length
		Extra_Vertices <- rep("NAME",No_Lost_Vertices)
		skip <- 0
		for(m in 1:length(unique(c(connections$V1, connections$V3)))){
			if(cluster_method!="highconn"){
				if((unique(c(connections$V1, connections$V3))[m] %in% V(Network)$name)==T){
				  skip <- skip + 1
				}else{
				  Extra_Vertices[m-skip] <- unique(c(connections$V1, connections$V3))[m]
				}
			}else{
				if((unique(c(connections$V1, connections$V3))[m] %in% gV)==T){
				  skip <- skip + 1
				}else{
				  Extra_Vertices[m-skip] <- unique(c(connections$V1, connections$V3))[m]
				}
			}
		}
	}else{
		Extra_Vertices <- NULL
	}

	if(Net_Length != 0){
		if(cluster_method=="markov"){
			Clustering <- mcl(x = as_adjacency_matrix(Network, type = "both", attr = "strength"), allow1 = T, addLoops = T, ESM = F, expansion = par_expansion, inflation = par_inflation, max.iter = 1000)
			if(is.atomic(Clustering)==T){
				print("An error occurred during the clustering")
			}else{
				No_Groups <- Clustering$K
			}
		}else if(cluster_method=="noclustering"){
			Clustering <- components(Network)
			No_Groups <- components(Network)$no
		}else if(cluster_method=="highconn"){
			Clustering <- highlyConnSG(Network)
			No_Groups <- length(Clustering$clusters)
		}else{
			if(cluster_method=="labelprop"){
				Clustering <- cluster_label_prop(Network, weights = E(Network)$strength)
			}else if(cluster_method=="community-edgebetweenness"){
				Clustering <- edge.betweenness.community(Network, weights = E(Network)$strength)
			}else if(cluster_method=="community-fastgreedy"){
				Clustering <- fastgreedy.community(Network, weights = E(Network)$strength)
			}else if(cluster_method=="community-spinglass"){
				Clustering <- spinglass.community(Network, weights = E(Network)$strength)
			}else if(cluster_method=="community-walktrap"){
				Clustering <- walktrap.community(Network, weights = E(Network)$strength)
			}
			No_Groups <- length(groups(Clustering))
		}

		Group_Members <- list()
		if(No_Groups == 0){
			Group_Members[[1]] <- NULL
		}else{
			for(j in 1:No_Groups){
				if(cluster_method=="markov"){
				  Group_Members[[j]] <- V(Network)$name[Clustering$Cluster==unique(Clustering$Cluster)[j]]
				}else if(cluster_method=="highconn"){
				  Group_Members[[j]] <- Clustering$clusters[[j]]
				}else{
				  Group_Members[[j]] <- groups(Clustering)[[j]]
				}
			}
		}
		Groups <- Group_Members
	}else{
		Groups <- list()
	}

	sink(args[3], append=T, split=F);
	AllGroups <- c(Groups,Extra_Vertices);
	for(G in 1:length(AllGroups)){
		for(S in 1:length(AllGroups[[G]])){
			cat(AllGroups[[G]][S],"\t");
		}
		cat("\n");
	}
	sink();
}
