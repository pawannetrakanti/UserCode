#! /bin/sh -f
echo "Running for dijet........... "

# Generate the response
#jet_response_analyzer_x  largeEta_binning_ptRebin.config -input dijet_pp_mergedpthatbins_Track8_Jet22MC.root -useweight true \
#jet_response_analyzer_x  fine_binning_ptRebin.config -input dijet_pp_mergedpthatbins_Track8_Jet22MC.root -useweight true \
#    -algs akVs3PFJetAnalyzer:0.3  akVs3CaloJetAnalyzer:0.3 ak3CaloJetAnalyzer:0.3 ak3PFJetAnalyzer:0.3 \
#    akVs4PFJetAnalyzer:0.4 akVs4CaloJetAnalyzer:0.4 ak4CaloJetAnalyzer:0.4 ak4PFJetAnalyzer:0.4 \
#    akVs5PFJetAnalyzer:0.5 akVs5CaloJetAnalyzer:0.5 ak5CaloJetAnalyzer:0.5 ak5PFJetAnalyzer:0.5 \
#    -output jra_hiF_ak_dijet.root
#

# Get the l3 corrections
jet_l3_correction_x -input jra_hiF_ak_dijet.root -era JEC_STARTHI53_LV1_Track8_Jet22_dijet -batch true -formats pdf \
    -algs \
    akVs3PFJetAnalyzer akVs4PFJetAnalyzer akVs5PFJetAnalyzer \
    ak3PFJetAnalyzer ak4PFJetAnalyzer ak5PFJetAnalyzer \
    akVs3CaloJetAnalyzer akVs4CaloJetAnalyzer akVs5CaloJetAnalyzer \
    ak3CaloJetAnalyzer ak4CaloJetAnalyzer ak5CaloJetAnalyzer \
    -output jra_hiF_ak_l3_dijet.root


# use the l3 output and response output to get l2 corrections
jet_l2_correction_x -input jra_hiF_ak_dijet.root -l3input jra_hiF_ak_l3_dijet.root -era JEC_STARTHI53_LV1_Track8_Jet22_dijet -batch true -formats pdf \
    -algs  akVs3CaloJetAnalyzer akVs4CaloJetAnalyzer akVs5CaloJetAnalyzer  \
    ak3CaloJetAnalyzer ak4CaloJetAnalyzer ak5CaloJetAnalyzer  \
    akVs3PFJetAnalyzer akVs4PFJetAnalyzer akVs5PFJetAnalyzer \
    ak3PFJetAnalyzer ak4PFJetAnalyzer ak5PFJetAnalyzer \
    -output jra_hiF_ak_l2_dijet.root
