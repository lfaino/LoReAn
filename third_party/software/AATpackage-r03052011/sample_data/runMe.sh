#!/bin/bash 

./cleanMe.pl

echo
echo
echo "*************************************************"
echo "running AAT nucleotide spliced alignment pipeline"
echo "*************************************************"
echo 
echo

cmd="../bin/AAT.pl -N -q arab.genomicSeq -s arab.cdna --dds '-f 100 -i 20 -o 75 -p 70 -a 2000' --filter '-c 10' --gap2 '-x 1' "

echo $cmd
eval $cmd


echo "**********************************************"
echo "running AAT protein spliced alignment pipeline"
echo "**********************************************"
echo
echo

 
cmd="../bin/AAT.pl -P -q arab.genomicSeq -s arab.pep --dps '-f 100 -i 30 -a 200' --filter '-c 10' --nap '-x 10' "
echo $cmd
eval $cmd

echo "************************************************"
echo "combining results into a multiple alignment file"
echo "************************************************"
echo
echo

cmd="../bin/show *gap2 *nap > multalignment.show.txt"
echo $cmd
eval $cmd

echo
echo
echo "done.  See 'multalignment.show.txt' to examine all spliced alignments."
echo



