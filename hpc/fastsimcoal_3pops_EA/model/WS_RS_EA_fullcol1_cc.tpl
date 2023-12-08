//Number of population samples (demes)
3 populations to simulate
//Population effective sizes (number of genes)
NWS$
NRS$
NEA$
//Sample sizes
30
16
10
//Growth rates
0
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 0 MIG02$
0 0 MIG12$
MIG20$ MIG21$ 0
//Migration matrix 1
0 MIG01C$ 0
MIG10C$ 0 0
0 0 0
//Migration matrix 2
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
5 historical event
T2 2 2 0 NEARES$ 0 1
T2 1 1 0 NRSRES$ 0 1
T2 0 0 0 NWSRES$ 0 1
T1 2 0 1 1 0 2
T1 1 0 1 RELANC$ 0 2
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 2.4e-9 OUTEXP