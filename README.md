Based on work from https://github.com/jiazhihao/HAG

build
make.sh

executables:
Run comparison of all three with ERGraphExperiment_optg2g1.sh \<p\> \<n\> where p=10*(ER graph edge probability) and n=ER graph node count. Parameter pair must exist in the synthetic graph sets folder. E.g., "ERGraphExperiment_optg2g1.sh 2 15" runs the experiment on ER graphs with 15 nodes, edge probability 0.2.
<p>
optimal solution:
bins.out
<p>
full greedy:
orig_singlelayer.out
<p>
partial greedy:
greedy2.out
greedy2order.out (ordered version does tie-breakers in a more intentional order to match up with original) 
<p>
hub heuristic:
bignode.out
<p>
pairs heuristic:
bignodepairs.out
