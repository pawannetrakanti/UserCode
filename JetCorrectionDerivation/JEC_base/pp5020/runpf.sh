#! /bin/sh -f
echo "Running for PF  dijet........... "

# Generate the response

jet_response_analyzer_x  jra_pf.config \
    -input /mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_20151229_ppReco_757p1_HcalRespCorrs_v4_00_mc.root \
    -useweight true \
    -output jra_hiF_akpf_ppReco_757p1_dijet.root \
    -algs \
    ak3PFJetAnalyzer:0.3 \
    ak4PFJetAnalyzer:0.3 \
    ak5PFJetAnalyzer:0.3 
exit
#    ak1PFJetAnalyzer:0.3 \
#    ak2PFJetAnalyzer:0.3 \
#    ak3PFJetAnalyzer:0.3 \
#    ak4PFJetAnalyzer:0.3 \
#    ak5PFJetAnalyzer:0.3 \
#    ak6PFJetAnalyzer:0.3 \
#    akPu1PFJetAnalyzer:0.3 \
#    akPu2PFJetAnalyzer:0.3 \
#    akPu3PFJetAnalyzer:0.3 \
#    akPu4PFJetAnalyzer:0.3 \
#    akPu5PFJetAnalyzer:0.3 \
#    akPu6PFJetAnalyzer:0.3 
#
#########

# Get the l3 corrections
jet_l3_correction_x -input jra_hiF_akpf_dijet.root -era JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3 -batch true -formats pdf \
   -output jra_hiF_akpf_l3_dijet.root \
   -algs \
   ak1PFJetAnalyzer ak2PFJetAnalyzer ak3PFJetAnalyzer ak4PFJetAnalyzer ak5PFJetAnalyzer ak6PFJetAnalyzer
#exit
##
# use the l3 output and response output to get l2 corrections
jet_l2_correction_x -input jra_hiF_akpf_dijet.root -l3input jra_hiF_akpf_l3_dijet.root -era JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3 -batch true -formats pdf \
    -output jra_hiF_akpf_l2_dijet.root \
    -algs \
    ak1PFJetAnalyzer ak2PFJetAnalyzer ak3PFJetAnalyzer ak4PFJetAnalyzer ak5PFJetAnalyzer  ak6PFJetAnalyzer
