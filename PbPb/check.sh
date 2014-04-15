# non Pu JEC on non Pu jets
#root -l -b <<EOF
#.L loadfwlite.C
#.x load_JEC.C 
#.x weight.C+("dijet_pp_mergedpthatbins_Track8_Jet22MC.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_Track8_Jet22MC_corrected_finelowpt.root",9999,1,true)
#.q
#EOF

# nonPu JEC on Pu jets in pp
#root -l -b <<EOF
#.L loadfwlite.C
#.x load_JEC.C 
#.x weight.C+("dijet_pp_mergedpthatbins_Track8_Jet22MC.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_Track8_Jet22MC_corrected_nonPuJets.root",9999,1,true)
#.q
#EOF

# Embedded pp
#root -l -b <<EOF
#.L loadfwlite.C
#.x load_JEC.C 
#.x weight.C+("dijet_pp_mergedpthatbins_embedpp.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_embedpp_corrected.root",9999,1,true)
#.q
#EOF


root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_check.C+("dijet_pp_mergedpthatbins_Track8_Jet22MC.root","/net/hidsk0001/d00/scratch/pawan/combinedPtHat/dijet_pp_mergedpthatbins_Track8_Jet22MC_corrected_finelowpt_onlymean.root",9999,1,true)
.q
EOF