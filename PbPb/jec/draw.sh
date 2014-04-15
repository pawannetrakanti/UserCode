#draw l2 corrections
#jet_draw_l2_correction_x -alg ak3PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg ak4PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg ak5PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg ak3CaloJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg ak4CaloJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg ak5CaloJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg akVs3PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg akVs4PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg akVs5PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg akVs3CaloJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg akVs4CaloJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#jet_draw_l2_correction_x -alg akVs5CaloJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles -outputFormat .pdf
#

jet_draw_l2_correction_x -alg ak3PFJetAnalyzer -filename jra_hiF_ak_l2_dijet.root -outputDir ./pdffiles/l2 -sameCanvas false -outputFormat .pdf

# draw l2l3 corrections 
#jet_draw_corrections_x -algs ak3pf -path "" -outputDir ./pdffiles/l2l3 -era JEC_STARTHI53_LV1_Track8_Jet22_dijet -useL2Cor true -useL3Cor true -outputFormat .pdf
