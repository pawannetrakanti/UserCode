root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt15/HiForest_v77_merged01/pt15_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp15.root",30,2.034E-01,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt30/HiForest_v77_merged01/pt30_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp30.root",50,1.075E-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt50/HiForest_v77_merged01/pt50_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp50.root",80,1.025E-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt80/HiForest_v77_merged01/pt80_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp80.root",120,9.865E-05,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt120/HiForest_v77_merged01/pt120_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp120.root",170,1.129E-05,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt170/HiForest_v77_merged01/pt170_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp170.root",220,1.465E-06,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt220/HiForest_v77_merged01/pt220_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp220.root",280,2.837E-07,false)
.q
EOF


root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_ppTracking.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt280/HiForest_v77_merged01/pt280_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp280.root",370,5.323E-08,false)
.q
EOF

list=`echo /net/hidsk0001/d00/scratch/pawan/combinedPtHat/ppTracking/pp*.root`
hadd /net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_prod16_ppTrackingMC.root $list
