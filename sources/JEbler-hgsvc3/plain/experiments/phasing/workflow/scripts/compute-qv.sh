#!/usr/bin/env bash

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./compute-qv.sh <read.meryl> <asm.fasta> <out> <threads> <memory>"
	echo
	echo -e "\t<read.meryl>:\tmeryl k-mer db of the (illumina) read set"
	echo -e "\t<asm.fasta>:\tassembly in FASTA format."
	echo -e "\t<out>:\tprefix of the output files."
	echo -e "\t<threads>:\tnumber of threads to use."
	echo -e "\t<memory>: max memory (GB) to be used for counting."
	echo
	echo "Computes the QV value of given assembly."
	exit 0
fi

read_db=$1
asm_fa=$2
name=$3
t=$4
m=$5


k=`meryl print $read_db | head -n 2 | tail -n 1 | awk '{print length($1)}'`
echo "Detected k-mer size $k"

echo "Count assembly k-mers ... "
meryl count k=$k threads=$t memory=$m output $name.meryl $asm_fa &> $name.meryl.log

echo "Determine assembly-specific k-mers ..."
meryl difference output $name.0.meryl $name.meryl $read_db &> $name.0.meryl.log

echo "Compute QV statistics ..."

ASM_ONLY=`meryl statistics $name.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
TOTAL=`meryl statistics $name.meryl  | head -n4 | tail -n1 | awk '{print $2}'`

if [[ $TOTAL -eq 0 ]]; then
  echo "[[ ERROR ]] :: assembly has no kmers."
  exit 1
else
  ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
  QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
  echo -e "$asm_fa\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" >> $name.qv
fi

rm -rf $name.meryl
rm -rf $name.0.meryl

echo "Done."
exit 0
