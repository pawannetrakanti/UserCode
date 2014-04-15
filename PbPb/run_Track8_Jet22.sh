root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT15_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp15.root",30,2.034E-01,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT30_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp30.root",50,1.075E-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT50_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp50.root",80,1.025E-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT80_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp80.root",120,9.865E-05,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT120_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp120.root",170,1.129E-05,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT170_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp170.root",220,1.465E-06,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT220_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp220.root",280,2.837E-07,false)
.q
EOF


root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT280_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp280.root",370,5.323E-08,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT370_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp370.root",460,5.934E-09,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT460_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp460.root",540,8.125E-10,false)
.q
EOF

root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight.C+("/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet22_cut1/hiForest_QCDpT540_STARTHI53_LV1_Track8_Jet22_1GeVcut.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp540.root",9999,1.467E-10,false)
.q
EOF

list=`echo /net/hidsk0001/d00/scratch/pawan/combinedPtHat/Track8_Jet22/pp*.root`
hadd /net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_Track8_Jet22MC.root $list
