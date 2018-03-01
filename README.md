# HitAnalyzer

ED analyzer to produce reco jet collection, pixel clusters and gen particles for b-tagging with nHits.
To run, change input file in config.py, then
Setup:
```
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
git clone git@github.com:thaarres/HitAnalyzer.git
cd HitAnalyzer
scram b
```
To run
```
cmsRun config.py
```
