#!/bin/bash

# Input
CLUSTER_FILE=$1
FCLUSTER_FILE=$2
TCG_ALIGNMENT=$3
SCORE_FILE=$4

# Output
NON_REC_PAIRS_CAT1=$5	#These pairs are from the same duplication and same fate
NON_REC_PAIRS_CAT2=$6	#These pairs are from different dups and different fates
REC_PAIRS=$7

CLUSTERS=`tail -n+2 $CLUSTER_FILE | cut -f2-`
IFS=$'\n'
CLUSTER_ONE_COUNTER=0
for CLUSTER_ONE in $CLUSTERS;
do
	#Check that this cluster is not broken up in the fate clustering
	PREV_SP_POS=0
	IFS=$'\t'
	for SPECIES in $CLUSTER_ONE;
	do
		NEXT_SP_POS=`grep -n $SPECIES $FCLUSTER_FILE | head -1 | cut -d ":" -f1`
		if [ "$NEXT_SP_POS" != "$PREV_SP_POS" ] && [ "$PREV_SP_POS" != "0" ]; then	#Then we will not use this cluster to make pairs
			continue 2
		fi
		PREV_SP_POS=$NEXT_SP_POS
	done

	SPECIES_ONE_COUNTER=0
	for SPECIES_ONE in $CLUSTER_ONE;
	do
		SPECIES_ONE_COUNTER=$((SPECIES_ONE_COUNTER+1))
		SPECIES_TWO_COUNTER=0
		for SPECIES_TWO in $CLUSTER_ONE;
		do
			SPECIES_TWO_COUNTER=$((SPECIES_TWO_COUNTER+1))
			if [ "$SPECIES_TWO_COUNTER" -gt "$SPECIES_ONE_COUNTER" ]; then #We only consider the upper triangle matrix
				#Get the gene names corresponding to species_one and species_two
				GENES=`egrep ">$SPECIES_ONE|>$SPECIES_TWO" $TCG_ALIGNMENT | cut -d " " -f1 | tr "\n" "\t"`
				TEST_SAME_FATE=`tail -n+2 $FCLUSTER_FILE | cut -f2- | egrep -c "$SPECIES_ONE|$SPECIES_TWO"`
				if [ "$TEST_SAME_FATE" == "2" ]; then
					echo -en ${CLUSTER_FILE::-7}"\t"$GENES"\n" >> $NON_REC_PAIRS_CAT1
				fi
			fi
		done
	done
	IFS=$'\n'

	CLUSTER_ONE_COUNTER=$((CLUSTER_ONE_COUNTER+1))
	CLUSTER_TWO_COUNTER=0
	for CLUSTER_TWO in $CLUSTERS;
	do
		#Check that this cluster is not broken up in the fate clustering
		PREV_SP_POS=0
		IFS=$'\t'
		for SPECIES in $CLUSTER_TWO;
		do
			NEXT_SP_POS=`grep -n $SPECIES $FCLUSTER_FILE | head -1 | cut -d ":" -f1`
			if [ "$NEXT_SP_POS" != "$PREV_SP_POS" ] && [ "$PREV_SP_POS" != "0" ]; then	#Then we will not use this cluster to make pairs
				continue 2
			fi
			PREV_SP_POS=$NEXT_SP_POS
		done

		CLUSTER_TWO_COUNTER=$((CLUSTER_TWO_COUNTER+1))
		if [ "$CLUSTER_TWO_COUNTER" -gt "$CLUSTER_ONE_COUNTER" ]; then #We only consider the upper triangle matrix
			#We now want all pairs between cluster_one and cluster_two
			IFS=$'\t'
			for SPECIES_ONE in $CLUSTER_ONE;
			do
				for SPECIES_TWO in $CLUSTER_TWO;
				do
					#Get the gene names corresponding to species_one and species_two
					GENES=`egrep ">$SPECIES_ONE|>$SPECIES_TWO" $TCG_ALIGNMENT | cut -d " " -f1 | tr "\n" "\t"`
					TEST_SAME_FATE=`tail -n+2 $FCLUSTER_FILE | cut -f2- | egrep -c "$SPECIES_ONE|$SPECIES_TWO"`
					if [ "$TEST_SAME_FATE" == "2" ]; then
						#If we only retrieve 2 fates we have to make sure that both species are present in the fcluster-file, otherwise these 2 fates apparently only include one of the species and we cannot say that the genes are recurrently differentiated.
						TEST_SPECIES_ONE=`tail -n+2 $FCLUSTER_FILE | cut -f2- | egrep -c "$SPECIES_ONE"`
						TEST_SPECIES_TWO=`tail -n+2 $FCLUSTER_FILE | cut -f2- | egrep -c "$SPECIES_TWO"`

						if [ "$TEST_SPECIES_ONE" == "2" ] && [ "$TEST_SPECIES_TWO" == "2" ]; then
							TEST_MAX_FATE_SP_ONE=`head -n-1 $SCORE_FILE | egrep -c "$SPECIES_ONE"`
							TEST_MAX_FATE_SP_TWO=`head -n-1 $SCORE_FILE | egrep -c "$SPECIES_TWO"`
							if [ "$TEST_MAX_FATE_SP_ONE" == "2" ] && [ "$TEST_MAX_FATE_SP_TWO" == "2" ]; then
								echo -en ${CLUSTER_FILE::-7}"\t"$GENES"\n" >> $REC_PAIRS
							fi
						else
							echo -e $CLUSTER_FILE":\tOne of the species is missing in the fclusters-file.\t"${GENES[*]}
						fi
					elif [ "$TEST_SAME_FATE" == "4" ]; then
						echo -en ${CLUSTER_FILE::-7}"\t"$GENES"\n" >> $NON_REC_PAIRS_CAT2
					else
						echo -e $CLUSTER_FILE":\tWarning: weird number of fates found in fclusters-file (perhaps both species are not present in it).\t"${GENES[*]}
					fi
				done
			done
			IFS=$'\n'
		fi
	done
done
