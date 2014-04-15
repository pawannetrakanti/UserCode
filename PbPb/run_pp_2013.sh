root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt15/HiForest_v81_merged01/pt15_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp15.root",30,2.034e-01,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt30/HiForest_v81_merged01/pt30_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp30.root",50,1.075e-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt50/HiForest_v81_merged01/pt50_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp50.root",80,1.025e-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt80/HiForest_v81_merged01/pt80_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp80.root",120,9.865e-05,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt120/HiForest_v81_merged01/pt120_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp120.root",170,1.129e-05,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt170/HiForest_v81_merged01/pt170_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp170.root",220,1.465e-06,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt220/HiForest_v81_merged01/pt220_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp220.root",280,2.837e-07,false)
.q
EOF


root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt280/HiForest_v81_merged01/pt280_pp2013_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp280.root",9999,5.323e-08,false)
.q
EOF

list=`echo /net/hidsk0001/d00/scratch/pawan/combinedPtHat/prod22v18/pp*.root`
hadd /net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_prod22v18MC.root $list
