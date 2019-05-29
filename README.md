#MyRandomWalks
These programs extends the TreeTraitMap class in the Beast-classic package (part of Beast2, 
a cross platform program for phylogenetic analysis).  The programs simulate various geographical 
population spreads including simple isotropic random walk, biased random walk, correlated random walk,
random walk with boundaries and long-range dispersal.  The program then performs Bayesian inference
using the MCMC sampling methods to estimate the location of the root of the phylogenetic tree from the 
locations of the leaves.

To run the programs dependencies have to be built between beast2, beast-classic, beast-geo, BeastLabs
and external library mason1.9

Each programme takes its corresponding examples/test...xml file as input argument and is run through the
main at beast.app.beastapp.BeastMain
