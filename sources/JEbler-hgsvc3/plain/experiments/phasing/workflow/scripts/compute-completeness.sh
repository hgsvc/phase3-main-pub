#!/usr/bin/env bash

if [ -z $1 ]; then
	echo "Usage: ./compute-completeness.sh <in.meryl> <asm.fasta> <out> <threads> <memory> "
	echo
	echo -e  "\t<in.meryl>:\tmeryl db of illumina data"
	echo -e "\t<asm.fasta>:\tassembly in FASTA format"
	echo -e "\t<out>:\tprefix of the output files"
	echo -e "\t<threads>:\t number of threads to use for counting"
	echo -e "\t<memory>:\t memory (GB) to use for counting"
	echo
	echo "Computes completeness of given assembly."
	exit -1
fi

db=$1
asm_fa=$2
name=$3
t=$4
m=$5


# Step 1 filter the read k-mers

k=`meryl print $db | head -n 2 | tail -n 1 | awk '{print length($1)}'`

echo "[compute-completeness] Generate histogram"
meryl histogram $db > $name.reads.hist

java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.jar $name.reads.hist > $name.reads.ploidy

cat $name.reads.ploidy

filt=`cat $name.reads.ploidy | grep -v "warning" | sed -n 2p | awk '{print $NF}'`

echo "[compute-completeness] Filter out kmers <= $filt"

meryl greater-than $filt output $name.reads.gt$filt.meryl $db
echo $filt > $db.filt


# Step 2 compute completeness of assembly
meryl count k=$k threads=$t memory=$m output $name.assembly.meryl $asm_fa &> $name.assembly.meryl.log
meryl intersect output $name.solid.meryl $name.assembly.meryl $name.reads.gt$filt.meryl

TOTAL=`meryl statistics $name.reads.gt$filt.meryl | head -n3 | tail -n1 | awk '{print $2}'`
ASM=`meryl statistics $name.solid.meryl | head -n3 | tail -n1 | awk '{print $2}'`
echo -e "${asm_fa}\tall\t${ASM}\t${TOTAL}" | awk '{print $0"\t"((100*$3)/$4)}' > $name.completeness.stats

rm -rf $name.solid.meryl
rm -rf $name.reads.gt$filt.meryl 
rm -rf $name.assembly.meryl

echo "[compute-completeness] Done"
exit 0
