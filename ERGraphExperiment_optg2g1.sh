p=$1
n=$2
kmax=3
testdir="syntheticGraphs/smallertest_n"$n"_p0p"$p"/"
outdir="ERGraphExperiment_n"$n"_p0p"$p"_k"$kmax"_optvsg1g2/"
readundir=1
mkdir $outdir
cumulativeoutfile=$outdir"actualcumulative.txt"
>$cumulativeoutfile
for id in {1..50}
do
	echo $id
	ingraph=$testdir"graph_"$id".txt"
	outfile=$outdir"actualscore"$id".txt"
	outfile2=$outdir"actualscorestep"$id".txt"
   tempout="tempout_bins"$id
   echo "(opt) ingraph: "$ingraph
#	>outfile
#	./bignodepairs.out $ingraph $kmax $readundir > $tempout
	./bins.out $ingraph $kmax $readundir  >$tempout
	cat $tempout | grep "cost saved by match:" |awk '{print $5}' > $outfile
	cat $tempout | grep "cost saved by match:" | awk '{print $5}' >> $cumulativeoutfile
	#cat $tempout | grep "stepscore" | awk '{print $2}' > $outfile2
#	cat $outfile | awk '{print $5}' >> $cumulativeoutfile
done


#cumulativeoutfile=$outdir"actualcumulative.txt"
cumulativeoutfileg2=$outdir"g2cumulative.txt"
>$cumulativeoutfileg2
for id in {1..50}
do
	echo $id
	ingraph=$testdir"graph_"$id".txt"
#	outfile=$outdir"actualscore"$id".txt"
#	outfile2=$outdir"actualscorestep"$id".txt"
    outfile=$outdir"g2score"$id".txt"
	outfile2=$outdir"g2scorestep"$id".txt"
   tempout="tempout_bins"$id
   echo "(g2) ingraph: "$ingraph
#	>outfile
#	./bignodepairs.out $ingraph $kmax $readundir > $tempout
	./greedy2.out $ingraph $kmax $readundir  >$tempout
	cat $tempout | grep "cost saved by match:" |awk '{print $5}' > $outfile
	cat $tempout | grep "cost saved by match:" | awk '{print $5}' >> $cumulativeoutfileg2
	#cat $tempout | grep "stepscore" | awk '{print $2}' > $outfile2
#	cat $outfile | awk '{print $5}' >> $cumulativeoutfile
done

cumulativeoutfileg1=$outdir"g1cumulative.txt"
>$cumulativeoutfileg1
for id in {1..50}
do
	echo $id
	ingraph=$testdir"graph_"$id".txt"
#	outfile=$outdir"actualscore"$id".txt"
#	outfile2=$outdir"actualscorestep"$id".txt"
    outfile=$outdir"g1score"$id".txt"
	outfile2=$outdir"g1scorestep"$id".txt"
   tempout="tempout_g1"$id
   echo "(g1) ingraph: "$ingraph
#	>outfile
#	./bignodepairs.out $ingraph $kmax $readundir > $tempout
	./orig_singlelayer.out $ingraph $kmax $readundir  >$tempout
	cat $tempout | grep "scoreSaved:" |awk '{print $2}' > $outfile
	cat $tempout | grep "scoreSaved:" | awk '{print $2}' >> $cumulativeoutfileg1
	#cat $tempout | grep "stepscore" | awk '{print $2}' > $outfile2
#	cat $outfile | awk '{print $5}' >> $cumulativeoutfile
done

echo "loops done"
twocol=$outdir"combinedcumulative.txt"
diffile1=$outdir"difs1.txt"
diffile2=$outdir"difs2.txt"
allres=$outdir"finalresults.txt"
paste $cumulativeoutfile $cumulativeoutfileg1 $cumulativeoutfileg2 > $twocol 
echo "paste done: "$cumulativeoutfile" "$cumulativeoutfileg1" "$cumulativeoutfileg2
cat $twocol | awk '{print $1-$2}' > $diffile1
cat $twocol |awk '{print $1-$3}' > $diffile2
paste $twocol $diffile1 $diffile2 > $allres
echo "opt < g1:"
cat $diffile1 | grep "\-" | wc -l
echo "opt = g1: "
cat $diffile1 | grep "0" | wc -l 
echo "opt < g2:"
cat $diffile2 | grep "\-" | wc -l
echo "opt = g2: "
cat $diffile2 | grep "0" | wc -l 




