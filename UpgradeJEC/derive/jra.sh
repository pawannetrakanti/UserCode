#!/bin/bash

doResponse=false
doCorrections=false
doApplyJEC=true

runalgo=ak4pf
#runalgo=ak4pfl2l3

if [ "$1" = "3" ]
then
##   Getting Jet response, resolution, fits etc...
    inputFile=JRA_ECAL-AGED-HCAL-AGED-SIPM-NOAGED_SCENARIO3.root
    outputFile=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3.root
    outputFitFile=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_f.root
    outputGraphFile=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_f_g.root

##   L3 corrections
#    inputForl3=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3.root
     inputForl3=$outputFitFile
     outputl3=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_l3.root

##   L2 corrections 
    inputl2=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_f.root
    inputl3=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_l3.root
    outputl2=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_l2.root

##  Check corrections
    jecpath=/home/pawan/t3store/UpgradeJEC/CMSSW_9_1_1_patch3/src/CondFormats/JetMETObjects/data

##  Apply JEC on fly
    outputFileJEC=Scene3/JRA_ECAL-AGED-HCAL-AGED-SIPM-NOAGED_SCENARIO3_JEC.root
#

## After jet energy corrections
    inputFileJEC=$outputFileJEC
    outputFile1=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_jec.root
    outputFitFileJEC=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_jec_f.root
    outputGraphFileJEC=Scene3/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_jec_f_g.root


elif [ "$1" = "4" ]
then
##   Getting Jet response, resolution, fits etc...
    inputFile=JRA_ECAL-AGED-HCAL-AGED-SIPM-AGED_SCENARIO4.root
    outputFile=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4.root
    outputFitFile=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4_f.root
    outputGraphFile=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4_f_g.root

##   L3 corrections    
    inputForl3=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4.root
    outputl3=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4_l3.root

##   L2 corrections     
    inputl2=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4.root
    inputl3=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4_l3.root
    outputl2=Scene4/jra_ecal-aged-hcal-aged-sipm-aged_scenario4_l2.root

##  Apply JEC on fly
    oututFileJEC=JRA_ECAL-AGED-HCAL-AGED-SIPM-AGED_SCENARIO4_JEC.root

else
##   Getting Jet response, resolution, fits etc...
    inputFile=JRA_ECAL-AGED-HCAL-NOAGED-SIPM-AGED_SCENARIO5.root
    outputFile=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5.root
    outputFitFile=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_f.root
    outputGraphFile=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_f_g.root

##   L3 corrections    
    inputForl3=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5.root
    outputl3=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_l3.root

##   L2 corrections    
    inputl2=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5.root
    inputl3=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_l3.root
    outputl2=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_l2.root

##  Apply JEC on fly
    outputFileJEC=Scene5/JRA_ECAL-AGED-HCAL-NOAGED-SIPM-AGED_SCENARIO5_JEC.root

# After jet energy corrections
    inputFileJEC=$outputFileJEC
    outputFile1=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_jec.root
    outputFitFileJEC=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_jec_f.root
    outputGraphFileJEC=Scene5/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_jec_f_g.root

fi

if [ "$doResponse" = true ]
then
# Generate the response
jet_response_analyzer_x $CMSSW_BASE/src/JetMETAnalysis/JetAnalyzers/config/jra_dr_coarsebinningeta_upgrade.config -input $inputFile -algs $runalgo:0.2 -output $outputFile -nbinsrelrsp 300 -relrspmax 6 -nbinsabsrsp 0 -nbinsetarsp 0 -nbinsphirsp 0 -etabarrelmin 0 -etabarrelmax 1.3 -etaiendcapmin 1.3 -etaiendcapmax 2.5 -etaoendcapmin 2.5 -etaoendcapmax 3.0 -etaforwardmin 3.0 -etaforwardmax 5.0 -doflavor false -flavorDefinition phys -useweight false

jet_response_fitter_x -input $outputFile 

jet_response_and_resolution_x -input $outputFitFile -dorelrsp true -algs $runalgo

# Graph resolution
#jet_inspect_graphs_x -inputs $outputGraphFile -algs $runalgo -variables RelResVsRefPt:JetEta@0
# Graph response
#jet_inspect_graphs_x -inputs $outputGraphFile -algs $runalgo -variables RelRspVsRefPt:JetEta@0 
fi

if [ "$doCorrections" = true ]
then
# Get the l3 corrections
    jet_l3_correction_x -input $inputForl3 -era PhaseIISummer16_25nsV3_MC_Scenario$1 -batch true -formats pdf \
	-output $outputl3 -algs $runalgo
    
# use the l3 output and response output to get l2 corrections
    jet_l2_correction_x -input $inputl2 -l3input $inputl3 -era PhaseIISummer16_25nsV3_MC_Scenario$1 -batch true -formats pdf \
	-output $outputl2 -algs $runalgo
fi


if [ "$doApplyJEC" = true ]
then
# Apply jec
#    jet_apply_jec_x -input $inputFile -era PhaseIISummer16_25nsV3_MC_Scenario$1  -levels 2 3 -output $outputFileJEC -algs $runalgo
    
#    runalgo=ak4pfl2l3
# After JEC is applied
#   jet_response_analyzer_x $CMSSW_BASE/src/JetMETAnalysis/JetAnalyzers/config/jra_dr_coarsebinningeta_upgrade.config -input $inputFileJEC -algs ${runalgo}l2l3:0.2 -output $outputFile1 -nbinsrelrsp 300 -relrspmax 6 -nbinsabsrsp 0 -nbinsetarsp 0 -nbinsphirsp 0 -etabarrelmin 0 -etabarrelmax 1.3 -etaiendcapmin 1.3 -etaiendcapmax 2.5 -etaoendcapmin 2.5 -etaoendcapmax 3.0 -etaforwardmin 3.0 -etaforwardmax 5.0 -doflavor false -flavorDefinition phys -useweight false
    
#   jet_response_fitter_x -input $outputFile1 -algs ${runalgo}l2l3
   
#   jet_response_and_resolution_x -input $outputFitFileJEC -dorelrsp true -algs ${runalgo}l2l3
    

# Graph resolution
#    jet_inspect_graphs_x -inputs $outputGraphFileJEC -algs ${runalgo}l2l3 -variables RelResVsRefPt:JetEta@0
#   jet_inspect_graphs_x -inputs $outputGraphFileJEC -algs $runalgo ${runalgo}l2l3 -variables RelResVsRefPt:JetEta@0 \
#       -colors 1 2 -markers 24 20
#       -removeFit true -xmin 35 -xmax 1500 -ymin 0.0 -ymax 0.4 \

# Graph response
    jet_inspect_graphs_x -inputs $outputGraphFileJEC -algs ${runalgo}l2l3 -variables RelRspVsRefPt:JetEta@1.0 
fi