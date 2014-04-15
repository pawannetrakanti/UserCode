#! /bin/sh -f
echo "Running jet response fitter for dijet........... "

# Jet response fitter
jet_response_fitter_x -input jra_hiF_ak_dijet.root -fittype 1 -verbose 1 -niter 10 \
    -algs \
    akVs3PFJetAnalyzer  
#    akVs4PFJetAnalyzer  akVs5PFJetAnalyzer 
#    akVs3CaloJetAnalyzer  akVs4CaloJetAnalyzer  akVs5CaloJetAnalyzer \
#    ak3PFJetAnalyzer  ak4PFJetAnalyzer  ak5PFJetAnalyzer \
#    ak3CaloJetAnalyzer  ak4CaloJetAnalyzer  ak5CaloJetAnalyzer 