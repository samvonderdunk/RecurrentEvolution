#!/bin/bash

# Input
FAMILY=$1		#Just give the name of the family, not any extension.
SCORE_FILE=$2
FCLUSTER_FILE=$3
PATTERN_FILE=$4

NON_REC_PAIRS_CAT1=$5	#These pairs are from the same duplication and same fate
NON_REC_PAIRS_CAT2=$6	#These pairs are from different dups and different fates
REC_PAIRS=$7

# Ouput
OUTPUT=$8

###
PAIR_CATEGORIES=($NON_REC_PAIRS_CAT1 $NON_REC_PAIRS_CAT2 $REC_PAIRS)
SCORE=`tail -1 $SCORE_FILE | cut -d " " -f3`
PAIRS=`grep -w $FAMILY $REC_PAIRS`

Z_scores=()
F_scores=()
PAIR_COUNT=0
SIG_PAIR_COUNT=0

if [ "$SCORE" == "1.0" ] || [ "$SCORE" == "0.0" ]; then
	echo -e "<ZF> = N/A" >> $OUTPUT
	exit
fi


IFS=$'\n'
for PAIR in $PAIRS;
do
	PAIR_COUNT=$((PAIR_COUNT+1))
	IFS=$' '
	GENES=(`echo $PAIR | cut -f2-`)

	#First check how the genes are ordered
	CHECK_FATE_G13=`egrep "${GENES[0]:1:10}|${GENES[2]:1:10}" $FCLUSTER_FILE | wc -l`
	CHECK_FATE_G14=`egrep "${GENES[0]:1:10}|${GENES[3]:1:10}" $FCLUSTER_FILE | wc -l`

	ORDER="OOOO"
	if [ "$CHECK_FATE_G13" == "1" ] && [ "$CHECK_FATE_G14" == "2" ]; then
		ORDER="ABAB"
	elif [ "$CHECK_FATE_G13" == "2" ] && [ "$CHECK_FATE_G14" == "1" ]; then
		ORDER="ABBA"
	fi

	PATTERNS=`cat $PATTERN_FILE | egrep "${GENES[0]:1:10}" | egrep "${GENES[2]:1:10}" | head -1 | cut -d " " -f5-`
	ABAB=`echo $PATTERNS | awk '{print $7+$19+$22+$35+$44+$45+$13+$31+$50+$8+$20}'`
	ABBA=`echo $PATTERNS | awk '{print $9+$17+$23+$30+$46+$47+$12+$36+$49+$10+$18}'`
	AABB=`echo $PATTERNS | awk '{print $4+$5+$14+$24+$26+$27+$29+$34+$41+$42+$51}'`
	if [ "$ORDER" == "ABBA" ]; then
		TEMP=$ABAB
		ABAB=$ABBA
		ABBA=$TEMP
	fi

	TOTAL=`echo $PATTERNS | awk '{print $1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30+$31+$32+$33+$34+$35+$36+$37+$38+$39+$40+$41+$42+$43+$44+$45+$46+$47+$48+$49+$50+$51+$52}'`
	dF=`echo $ABAB $ABBA | awk '{print $1-$2}'`
	SEMI_TOTAL=`echo $ABAB $ABBA $AABB | awk '{print $1+$2+$3}'`

	if [ "$SEMI_TOTAL" == "0" ]; then
		F=0
		VdF1=0
		ZF1=0
	else
		F=`echo $ABAB $ABBA $SEMI_TOTAL | awk '{print ($1-$2)/$3}'`
		VdF1=`echo $ABAB $ABBA $SEMI_TOTAL | awk '{print $3*($1/$3)*(1-($1/$3)) + $3*($2/$3)*(1-($2/$3)) + 2*$3*($1/$3)*($2/$3)}'`
		if [ "$VdF1" == "0" ]; then
			ZF1=0
		else
			ZF1=`echo $dF $VdF1 | awk '{print $1/(sqrt($2))}'`
		fi
	fi

	VdF2=`echo $PATTERNS $TOTAL | awk '{print $7" "$19" "$22" "$35" "$44" "$45" "$13" "$31" "$50" "$8" "$20" "$9" "$17" "$23" "$30" "$46" "$47" "$12" "$36" "$49" "$10" "$18" "$53}' | awk '{Vsum=0; for (i=1; i<=22; i++) Vsum+= $23*($i/$23)*(1-($i/$23)); for (i=1; i<=22; i++) for (j=i+1; j<=22; j++) Vsum+= 2*$23*($i/$23)*($j/$23); print Vsum}'`
	if [ "$VdF2" == "0" ]; then
		continue	#Choose to not include this pair, because its Z-score is basically not defined.
	else
		ZF2=`echo $dF $VdF2 | awk '{print $1/(sqrt($2))}'`
	fi

	#ZF2 seems only slightly lower than ZF1 for each pair, so let's use the slightly more conservative ZF2.
	F_scores+=($F)
	Z_scores+=($ZF2)

	if (( $(echo "$ZF2 > 1.96" |bc -l) )); then
		SIG_PAIR_COUNT=$((SIG_PAIR_COUNT+1))
	fi
	
done

IFS=$' '
MIN_Z=100
MAX_Z=0
for Z in ${Z_scores[*]}; do
	if (( $(echo "$Z > $MAX_Z" |bc -l) )); then
		MAX_Z=$Z
	fi

	if (( $(echo "$Z < $MIN_Z" |bc -l) )); then
		MIN_Z=$Z
	fi
done

MIN_F=100
MAX_F=0
for F in ${F_scores[*]}; do
	if (( $(echo "$F > $MAX_F" |bc -l) )); then
		MAX_F=$F
	fi

	if (( $(echo "$F < $MIN_F" |bc -l) )); then
		MIN_F=$F
	fi
done

Z_LENGTH=${#Z_scores[@]}
if ! (($Z_LENGTH%2)); then
	HALF_Z_LENGTH=$((Z_LENGTH/2))
	HALF_Z_LENGTH=$((HALF_Z_LENGTH+1))
	MEDIAN_Z=`echo ${Z_scores[*]} | tr " " "\n" | sort -n | head -$HALF_Z_LENGTH | tail -2 | awk '{sum+=$1} END {print sum/2}'`
	MEDIAN_F=`echo ${F_scores[*]} | tr " " "\n" | sort -n | head -$HALF_Z_LENGTH | tail -2 | awk '{sum+=$1} END {print sum/2}'`
else
	HALF_Z_LENGTH=$((Z_LENGTH/2))
	HALF_Z_LENGTH=$((HALF_Z_LENGTH+1))
	MEDIAN_Z=`echo ${Z_scores[*]} | tr " " "\n" | sort -n | head -$HALF_Z_LENGTH | tail -1`
	MEDIAN_F=`echo ${F_scores[*]} | tr " " "\n" | sort -n | head -$HALF_Z_LENGTH | tail -1`
fi

echo $PAIR_COUNT $SIG_PAIR_COUNT $MIN_Z $MEDIAN_Z $MAX_Z ${Z_scores[*]} | awk '{sum=0; for (i=6; i<=NF; i++) sum+=$i; printf ("<ZF> = %.3f\t(%i/%i significant pairs)\n", sum/(NF-5), $2, $1)}' >> $OUTPUT
