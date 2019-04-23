# HepSim/auroratest

This is benchmark example showing  IO and CPU intensive 
calcutaions which are typical for collider experiments.
This program runs over ProMC files on a Linux/bash system using the [HepSim Monte Carlo repository](http://atlaswww.hep.anl.gov/hepsim/), and then 
builds anti-KT jets (exclusive 6-jet definition) and k-means 
clusters (6 clusters) in the eta-phi space.
The k-means algorithm runs in the Euclidian space using 300 passes. 
Note that the centers of anti-kT jets and the positions of the k-means clusters are not expected to be
the same. 

 1. Install ProMC (http://atlaswww.hep.anl.gov/asc/promc/) and ROOT
 2. Check the installation. The variables: 

```
   echo $PROMC
   echo $ROOTSYS
```
  they should return the installation paths. 

 3. Compile as "make"

 4. Download ProMC files from HepSim and put them to the "data" directory. Use hs-tools as: 
  
``` 
   hs-get tev100_higgs_ttbar_mg5 data
```
   See the HepSim documentation. 

 5. Process all files inside the directory "data" using the command "./example".

 6. Loook at the output root file "output.root" with histograms.


References:

 1. FastJet (fjcore v3.3.2). M. Cacciari, G.P. Salam and G. Soyez, Eur.Phys.J. C72 (2012) 1896 [arXiv:1111.6097]  

 2. k-means for HEP studies: Chekanov, S. Eur. Phys. J  C 47 (2006), p. 611. ANL-HEP-PR-05-118   E-print: hep-ph/0512027

 3. ProIO: S.Chekanov, E.May, K. Strand, P. Van Gemmeren, ProMC: Input-output data format for HEP applications using varint encoding, ANL-HEP-PR-13-41, arXiv:1311.1229, Computer Physics Communications 185 (2014), pp. 2629-2635


S.Chekanov (ANL) 
