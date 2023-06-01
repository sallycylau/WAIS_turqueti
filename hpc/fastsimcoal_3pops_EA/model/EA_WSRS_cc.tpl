//Number of population samples (demes)
3 populations to simulate
//Population effective sizes (number of genes)
NEA
NWS
NRS
//Sample sizes
10
30
16
//Growth rates
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG01$ MIG02$
MIG10$ 0 0
MIG20$ 0 0
//Migration matrix 1
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
2 historical event
T2 2 1 1 1 0 1
T1 1 0 1 RELANC 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 2.4e-9 OUTEXP
