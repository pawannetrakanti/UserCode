root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_embedpp.C+("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/embedpp/pp30.root",50,1.075e-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_embedpp.C+("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/embedpp/pp50.root",80,1.025e-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_embedpp.C+("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/embedpp/pp80.root",100,9.865e-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_embedpp.C+("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0/0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/embedpp/pp100.root",120,3.069e-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_embedpp.C+("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0/0.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/embedpp/pp120.root",9999,1.129E-05,false)
.q
EOF

list=`echo /net/hidsk0001/d00/scratch/pawan/combinedPtHat/embedpp/pp*.root`
hadd /net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_embedpp.root $list
