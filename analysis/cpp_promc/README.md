# HepSim/analysis/cpp-promc

This example shows how to process ProMC files on a Linux/bash system using the [HepSim Monte Carlo repository](http://atlaswww.hep.anl.gov/hepsim/). The example builds anti-KT jets, and fill ROOT histograms.


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

S.Chekanov (ANL) 
