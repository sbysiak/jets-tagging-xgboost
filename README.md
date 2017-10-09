# Jets tagging with XGBoost
---------------------------

Analysis of possibility of using boosted decision trees (XGBoost implementation) for tagging jets flavour,  
done for my Bachelor thesis. 

For now it is possible to generate events with Monte Carlo generator [Pythia8](http://home.thep.lu.se/~torbjorn/Pythia.html) and reconstruct them with [FastJet](http://fastjet.fr/) (by default with _anti-kt_ algorithm _R_=0.4)

## Prerequisites
In order to run any of generation & reconstruction parts, having ALICE software (AliROOT, AliPhysics) installed is required.  
It is accessible via aliBuild-tool:

[https://dberzano.github.io/alice/alibuild/](https://dberzano.github.io/alice/alibuild/ "Dario Berzano: Getting started with aliBuild")

In case you installed required software in any other way, note having all environment variables ($ALICE_ROOT, FASTJET_ROOT ...) set correctly and that the following structure of alice directories is asserted:
```
-- $HOME
    |-- alice
        |-- ali-master
        |-- sw
```
otherwise, a few changes in paths included in in run_\*.C will be needed.


Whole repository was developed with following versions of packages:
```
  1) BASE/1.0
  2) AliEn-Runtime/v2-19-le-1
  3) GSL/v1.16-1
  4) ROOT/v5-34-30-alice-1
  5) boost/v1.59.0-1
  6) cgal/v4.6.3-1
  7) fastjet/v3.1.3_1.020-1
  8) GEANT3/v2-1-1
  9) GEANT4/v4.10.01.p03-1
 10) GEANT4_VMC/v3-2-p1-1
 11) AliRoot/0-1
 12) AliPhysics/latest-ali-master
 ```
  
  

## Getting started
Any part of code written in (Ali)ROOT 
(basically concerning generation and reconstruction) 
should be run inside _alienv_:
```
cd $HOME/alice/sw/; 
sudo alienv enter AliPhysics/latest-ali-master -a <ARCHITECTURE, e.g. ubuntu1404_x86-64>;
```
### Event generation & reconstruction
In order to generate first events using [Pythia8](http://home.thep.lu.se/~torbjorn/Pythia.html) and reconstruct them with [FastJet](http://fastjet.fr/), type:
```
aliroot 'generation/run_gener_recon.C(100,2,50,5)'
```


It is also possible to separately execute generation and reconstruction steps by:
```
aliroot 'generation/run_gener.C(100,2,50,5)'
```
and
```
aliroot 'generation/run_recon.C(100,2,50,5)'
```

### Parameters
For each __run\_\*.C__ macro the same set of parameters is used:
```
aliroot 'generation/run_gener_recon.C(nev, type, parton_en, q_id, [longFileName])'
```
  * __nev__           -- number of events
  * __type__          -- 1 for gluon-gluon and 2 for quark-antiquark pair
  * __parton_en__     -- energy of quark or gluon in GeV
  * __q_id__ 		  -- quark flavour: 1=u, 2=d, 3=s, 4=c, 5=b, 6=t;
        ignored for   type=1
  * __longFileName__  -- [optional, default: _"kTRUE"_] include params above in output files names;
  if _longFileName_ == _"kTRUE"_ then output files will contain parameters, like _Kinematics_g_en50_nev100.root_
  
For example:
```
aliroot 'generation/run_gener_recon.C(500,2,50,5,"kFALSE")'
```
runs both generation and reconstruction of __500__ events, each consisting of 50 GeV ___b -- anti-b___ __quarks__ pair.
while
```
aliroot 'generation/run_gener.C(1000,1,150,0,"kTRUE")'
```
will generate __1000__ events of __gluon--gluon__ pairs with _E/gluon_ = 150 GeV (without reconstruction).

 ### Outputs
 Macros generating events produce two files: _Kinematics.root_ - contaning event-by-event kinematics information about all particles produced and _galice.root_ - contaning meta - information about events.
 
 Macros reconstructing events produce: _histos.root_ - with some auxiliary histograms and  ___treeOfJets.root___ - in fact __most important__ file, which contains ROOT tree with four branches: 
 * jetsBranch -- information concerning jets originated directly from _FastJet_
 * evParticlesBranch -- information about all particles processed in the events
 * jetParticlesBranch -- information about particles which construct jets (according to certain reconstruction algorithm)
 * jetObservablesBranch -- high-level observables calculated for jets:
   * ptRel - relative pt of the electron
   * eIn - no. electrons in jets (usually binary information if _e-_ is present in jet)
   * mult - multiplicity
   * radMom - radial moment
   * angular - angularity
   * svR - distance between primary and secondary vertices
   * svN - no. particles originated from secondary vertex


## Performance
Execution time of generation and reconstruction increases 
improportionally to number of events. Couple examples of run times 
on average laptop:
  * no. events   : time
  - 300 events   : 16 sec
  - 1000 events  : 31 sec
  - 3000 events  : 1 min 30 sec
  - 10000 events : 7 min 40 sec
  - 30000 events : 54 min

In order to simulate larger numbers of events, consider merging with [hadd](https://root.cern.ch/root/HowtoMerge.html "ROOT: hadd") a couple of files, each containing fraction of desired no. events.
