#! /bin/sh 
echo "Running for Calo dijet........... "

# Generate the response

jet_response_analyzer_x  jra_calo.config -input HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_20160104_HIReco_757p1_HcalRespCorrs_v4_00_mc.root -useweight true \
    -output jra_hiF_akcalo_HIReco_757p1_dijet.root \
    -algs \
    ak1CaloJetAnalyzer:0.3 \
    ak2CaloJetAnalyzer:0.3 \
    ak3CaloJetAnalyzer:0.3 \
    ak4CaloJetAnalyzer:0.3 \
    ak5CaloJetAnalyzer:0.3 \
    ak6CaloJetAnalyzer:0.3 \
    akPu1CaloJetAnalyzer:0.3 \
    akPu2CaloJetAnalyzer:0.3 \
    akPu3CaloJetAnalyzer:0.3 \
    akPu4CaloJetAnalyzer:0.3 \
    akPu5CaloJetAnalyzer:0.3 \
    akPu6CaloJetAnalyzer:0.3 
exit
###########

# Get the l3 corrections
jet_l3_correction_x -input jra_hiF_akcalo_dijet.root -era JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3 -batch true -formats pdf \
   -output jra_hiF_akcalo_l3_dijet.root \
   -algs \
   ak1CaloJetAnalyzer ak2CaloJetAnalyzer ak3CaloJetAnalyzer ak4CaloJetAnalyzer ak5CaloJetAnalyzer ak6CaloJetAnalyzer 
#exit
###
# use the l3 output and response output to get l2 corrections
jet_l2_correction_x -input jra_hiF_akcalo_dijet.root -l3input jra_hiF_akcalo_l3_dijet.root -era JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3 -batch true -formats pdf \
    -output jra_hiF_akcalo_l2_dijet.root \
    -algs \
    ak1CaloJetAnalyzer ak2CaloJetAnalyzer ak3CaloJetAnalyzer ak4CaloJetAnalyzer ak5CaloJetAnalyzer  ak6CaloJetAnalyzer 
