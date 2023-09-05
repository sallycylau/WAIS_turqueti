//Number of population samples (demes)
4 populations to simulate
//Population effective sizes (number of genes)
NWS$
NAS$
NRS$
NEA$
//Sample sizes
30
8
16
10
//Growth rates
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG01$ 0 MIG03$
MIG10$ 0 MIG12$ 0
0 MIG21$ 0 MIG23$
MIG30$ 0 MIG32$ 0
//Migration matrix 1
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
9 historical event
T3 3 3 0 NEARES$ 0 1
T3 2 2 0 NRSRES$ 0 1
T3 1 1 0 NASRES$ 0 1
T3 0 0 0 NWSRES$ 0 1
T2 3 0 1 1 0 1
T2 2 0 1 1 0 1
T2 1 0 1 RES1$ 0 1
T1 0 0 0 RELANC$ 0 1
T0 0 0 0 RELANC2$ 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 2.4e-9 OUTEXP