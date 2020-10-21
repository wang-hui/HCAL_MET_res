# MET response and resolution
This is a git repo for MET study using nanoAOD ntuples  

1. setup CMSSW. Any 10xx should work
```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
```

2. checkout nanoAOD-tools
```
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b -j 8
```

3. checkout this repo
```
git clone https://github.com/wang-hui/HCAL_MET_res.git
```

4. local test
```
cd HCAL_MET_res
python HCAL_MET_res.py FileList/DoubleMuon_Run2018A_NANOAOD.list
```
