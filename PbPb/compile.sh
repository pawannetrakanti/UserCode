root -l -b<<EOF
.x loadfwlite.C
.L FWCore/Utilities/src/typelookup.cc++
.L CondFormats/JetMETObjects/src/Utilities.cc++
.L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc++
.L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc++
.L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc++
.L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc++
.L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc++
EOF
