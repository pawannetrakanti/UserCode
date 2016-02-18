//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 18 07:26:09 2015 by ROOT version 6.05/03
// from TTree HltTree/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODpp5TeV/HiForestAOD_withTupel_pp_MC_Z30mumuJet_v1.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

class Hlts {
public :
   Hlts(){};
   ~Hlts(){};

   // Declaration of leaf types
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[10];   //[NL1IsolEm]
   Float_t         L1IsolEmE[10];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[10];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[10];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[10];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[10];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[10];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[10];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[10];   //[NL1Mu]
   Float_t         L1MuE[10];   //[NL1Mu]
   Float_t         L1MuEta[10];   //[NL1Mu]
   Float_t         L1MuPhi[10];   //[NL1Mu]
   Int_t           L1MuIsol[10];   //[NL1Mu]
   Int_t           L1MuMip[10];   //[NL1Mu]
   Int_t           L1MuFor[10];   //[NL1Mu]
   Int_t           L1MuRPC[10];   //[NL1Mu]
   Int_t           L1MuQal[10];   //[NL1Mu]
   Int_t           L1MuChg[10];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[10];   //[NL1CenJet]
   Float_t         L1CenJetE[10];   //[NL1CenJet]
   Float_t         L1CenJetEta[10];   //[NL1CenJet]
   Float_t         L1CenJetPhi[10];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[10];   //[NL1ForJet]
   Float_t         L1ForJetE[10];   //[NL1ForJet]
   Float_t         L1ForJetEta[10];   //[NL1ForJet]
   Float_t         L1ForJetPhi[10];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[10];   //[NL1Tau]
   Float_t         L1TauE[10];   //[NL1Tau]
   Float_t         L1TauEta[10];   //[NL1Tau]
   Float_t         L1TauPhi[10];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1EtTot;
   Float_t         L1Mht;
   Float_t         L1MhtPhi;
   Float_t         L1EtHad;
   Int_t           L1HfRing1EtSumPositiveEta;
   Int_t           L1HfRing2EtSumPositiveEta;
   Int_t           L1HfRing1EtSumNegativeEta;
   Int_t           L1HfRing2EtSumNegativeEta;
   Int_t           L1HfTowerCountPositiveEtaRing1;
   Int_t           L1HfTowerCountNegativeEtaRing1;
   Int_t           L1HfTowerCountPositiveEtaRing2;
   Int_t           L1HfTowerCountNegativeEtaRing2;
   ULong64_t       Event;
   Int_t           LumiBlock;
   Int_t           Run;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           digitisation_step;
   Int_t           digitisation_step_Prescl;
   Int_t           L1simulation_step;
   Int_t           L1simulation_step_Prescl;
   Int_t           digi2raw_step;
   Int_t           digi2raw_step_Prescl;
   Int_t           HLTriggerFirstPath;
   Int_t           HLTriggerFirstPath_Prescl;
   Int_t           HLT_Physics_v1;
   Int_t           HLT_Physics_v1_Prescl;
   Int_t           HLT_Random_v1;
   Int_t           HLT_Random_v1_Prescl;
   Int_t           HLT_ZeroBias_v1;
   Int_t           HLT_ZeroBias_v1_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           raw2digi_step;
   Int_t           raw2digi_step_Prescl;
   Int_t           L1Reco_step;
   Int_t           L1Reco_step_Prescl;
   Int_t           L1_AlwaysTrue;
   Int_t           L1_AlwaysTrue_Prescl;
   Int_t           L1_DoubleEG6_HTT150;
   Int_t           L1_DoubleEG6_HTT150_Prescl;
   Int_t           L1_DoubleEG_15_10;
   Int_t           L1_DoubleEG_15_10_Prescl;
   Int_t           L1_DoubleEG_22_10;
   Int_t           L1_DoubleEG_22_10_Prescl;
   Int_t           L1_DoubleIsoTau28er;
   Int_t           L1_DoubleIsoTau28er_Prescl;
   Int_t           L1_DoubleIsoTau32er;
   Int_t           L1_DoubleIsoTau32er_Prescl;
   Int_t           L1_DoubleIsoTau36er;
   Int_t           L1_DoubleIsoTau36er_Prescl;
   Int_t           L1_DoubleIsoTau40er;
   Int_t           L1_DoubleIsoTau40er_Prescl;
   Int_t           L1_DoubleMu0;
   Int_t           L1_DoubleMu0_Prescl;
   Int_t           L1_DoubleMu0_Eta1p6_WdEta18;
   Int_t           L1_DoubleMu0_Eta1p6_WdEta18_Prescl;
   Int_t           L1_DoubleMu0_Eta1p6_WdEta18_OS;
   Int_t           L1_DoubleMu0_Eta1p6_WdEta18_OS_Prescl;
   Int_t           L1_DoubleMu6_EG6;
   Int_t           L1_DoubleMu6_EG6_Prescl;
   Int_t           L1_DoubleMu7_EG7;
   Int_t           L1_DoubleMu7_EG7_Prescl;
   Int_t           L1_DoubleMuOpen;
   Int_t           L1_DoubleMuOpen_Prescl;
   Int_t           L1_DoubleMu_10_0_WdEta18;
   Int_t           L1_DoubleMu_10_0_WdEta18_Prescl;
   Int_t           L1_DoubleMu_10_3p5;
   Int_t           L1_DoubleMu_10_3p5_Prescl;
   Int_t           L1_DoubleMu_10_Open;
   Int_t           L1_DoubleMu_10_Open_Prescl;
   Int_t           L1_DoubleMu_12_5;
   Int_t           L1_DoubleMu_12_5_Prescl;
   Int_t           L1_DoubleJetC100;
   Int_t           L1_DoubleJetC100_Prescl;
   Int_t           L1_DoubleJetC112;
   Int_t           L1_DoubleJetC112_Prescl;
   Int_t           L1_DoubleJetC52;
   Int_t           L1_DoubleJetC52_Prescl;
   Int_t           L1_DoubleJetC56_ETM60;
   Int_t           L1_DoubleJetC56_ETM60_Prescl;
   Int_t           L1_DoubleJetC84;
   Int_t           L1_DoubleJetC84_Prescl;
   Int_t           L1_DoubleTauJet40er;
   Int_t           L1_DoubleTauJet40er_Prescl;
   Int_t           L1_EG25er_HTT100;
   Int_t           L1_EG25er_HTT100_Prescl;
   Int_t           L1_ETM100;
   Int_t           L1_ETM100_Prescl;
   Int_t           L1_ETM30;
   Int_t           L1_ETM30_Prescl;
   Int_t           L1_ETM40;
   Int_t           L1_ETM40_Prescl;
   Int_t           L1_ETM50;
   Int_t           L1_ETM50_Prescl;
   Int_t           L1_ETM60;
   Int_t           L1_ETM60_Prescl;
   Int_t           L1_ETM70;
   Int_t           L1_ETM70_Prescl;
   Int_t           L1_ETT15_BptxAND;
   Int_t           L1_ETT15_BptxAND_Prescl;
   Int_t           L1_ETT40;
   Int_t           L1_ETT40_Prescl;
   Int_t           L1_ETT60;
   Int_t           L1_ETT60_Prescl;
   Int_t           L1_HTT100;
   Int_t           L1_HTT100_Prescl;
   Int_t           L1_HTT125;
   Int_t           L1_HTT125_Prescl;
   Int_t           L1_HTT150;
   Int_t           L1_HTT150_Prescl;
   Int_t           L1_HTT175;
   Int_t           L1_HTT175_Prescl;
   Int_t           L1_HTT200;
   Int_t           L1_HTT200_Prescl;
   Int_t           L1_HTT250;
   Int_t           L1_HTT250_Prescl;
   Int_t           L1_HTT75;
   Int_t           L1_HTT75_Prescl;
   Int_t           L1_IsoEG20er_TauJet20er_NotWdEta0;
   Int_t           L1_IsoEG20er_TauJet20er_NotWdEta0_Prescl;
   Int_t           L1_Mu0er_ETM40;
   Int_t           L1_Mu0er_ETM40_Prescl;
   Int_t           L1_Mu0er_ETM55;
   Int_t           L1_Mu0er_ETM55_Prescl;
   Int_t           L1_Mu10er_ETM30;
   Int_t           L1_Mu10er_ETM30_Prescl;
   Int_t           L1_Mu10er_ETM50;
   Int_t           L1_Mu10er_ETM50_Prescl;
   Int_t           L1_Mu12_EG10;
   Int_t           L1_Mu12_EG10_Prescl;
   Int_t           L1_Mu14er_ETM30;
   Int_t           L1_Mu14er_ETM30_Prescl;
   Int_t           L1_Mu16er_IsoTau28er;
   Int_t           L1_Mu16er_IsoTau28er_Prescl;
   Int_t           L1_Mu16er_IsoTau32er;
   Int_t           L1_Mu16er_IsoTau32er_Prescl;
   Int_t           L1_Mu16er_TauJet20er;
   Int_t           L1_Mu16er_TauJet20er_Prescl;
   Int_t           L1_Mu20_EG10;
   Int_t           L1_Mu20_EG10_Prescl;
   Int_t           L1_Mu3_JetC16_WdEtaPhi2;
   Int_t           L1_Mu3_JetC16_WdEtaPhi2_Prescl;
   Int_t           L1_Mu3_JetC52_WdEtaPhi2;
   Int_t           L1_Mu3_JetC52_WdEtaPhi2_Prescl;
   Int_t           L1_Mu4_EG18;
   Int_t           L1_Mu4_EG18_Prescl;
   Int_t           L1_Mu5_DoubleEG5;
   Int_t           L1_Mu5_DoubleEG5_Prescl;
   Int_t           L1_Mu5_EG15;
   Int_t           L1_Mu5_EG15_Prescl;
   Int_t           L1_Mu5_EG20;
   Int_t           L1_Mu5_EG20_Prescl;
   Int_t           L1_Mu5_IsoEG18;
   Int_t           L1_Mu5_IsoEG18_Prescl;
   Int_t           L1_Mu6_DoubleEG10;
   Int_t           L1_Mu6_DoubleEG10_Prescl;
   Int_t           L1_Mu6_HTT100;
   Int_t           L1_Mu6_HTT100_Prescl;
   Int_t           L1_Mu8_HTT50;
   Int_t           L1_Mu8_HTT50_Prescl;
   Int_t           L1_QuadMu0;
   Int_t           L1_QuadMu0_Prescl;
   Int_t           L1_QuadJetC36_TauJet52;
   Int_t           L1_QuadJetC36_TauJet52_Prescl;
   Int_t           L1_QuadJetC40;
   Int_t           L1_QuadJetC40_Prescl;
   Int_t           L1_QuadJetC60;
   Int_t           L1_QuadJetC60_Prescl;
   Int_t           L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1;
   Int_t           L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1_Prescl;
   Int_t           L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1;
   Int_t           L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1_Prescl;
   Int_t           L1_SingleEG10;
   Int_t           L1_SingleEG10_Prescl;
   Int_t           L1_SingleEG15;
   Int_t           L1_SingleEG15_Prescl;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG25;
   Int_t           L1_SingleEG25_Prescl;
   Int_t           L1_SingleEG2_BptxAND;
   Int_t           L1_SingleEG2_BptxAND_Prescl;
   Int_t           L1_SingleEG30;
   Int_t           L1_SingleEG30_Prescl;
   Int_t           L1_SingleEG35;
   Int_t           L1_SingleEG35_Prescl;
   Int_t           L1_SingleEG40;
   Int_t           L1_SingleEG40_Prescl;
   Int_t           L1_SingleEG5;
   Int_t           L1_SingleEG5_Prescl;
   Int_t           L1_SingleIsoEG18er;
   Int_t           L1_SingleIsoEG18er_Prescl;
   Int_t           L1_SingleIsoEG20;
   Int_t           L1_SingleIsoEG20_Prescl;
   Int_t           L1_SingleIsoEG20er;
   Int_t           L1_SingleIsoEG20er_Prescl;
   Int_t           L1_SingleIsoEG22er;
   Int_t           L1_SingleIsoEG22er_Prescl;
   Int_t           L1_SingleIsoEG25;
   Int_t           L1_SingleIsoEG25_Prescl;
   Int_t           L1_SingleIsoEG25er;
   Int_t           L1_SingleIsoEG25er_Prescl;
   Int_t           L1_SingleIsoEG30er;
   Int_t           L1_SingleIsoEG30er_Prescl;
   Int_t           L1_SingleMu12;
   Int_t           L1_SingleMu12_Prescl;
   Int_t           L1_SingleMu14er;
   Int_t           L1_SingleMu14er_Prescl;
   Int_t           L1_SingleMu16;
   Int_t           L1_SingleMu16_Prescl;
   Int_t           L1_SingleMu16er;
   Int_t           L1_SingleMu16er_Prescl;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu20_Prescl;
   Int_t           L1_SingleMu20er;
   Int_t           L1_SingleMu20er_Prescl;
   Int_t           L1_SingleMu25;
   Int_t           L1_SingleMu25_Prescl;
   Int_t           L1_SingleMu30;
   Int_t           L1_SingleMu30_Prescl;
   Int_t           L1_SingleMu5;
   Int_t           L1_SingleMu5_Prescl;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu7_Prescl;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuOpen;
   Int_t           L1_SingleMuOpen_Prescl;
   Int_t           L1_SingleMuOpen_NotBptxOR;
   Int_t           L1_SingleMuOpen_NotBptxOR_Prescl;
   Int_t           L1_SingleJet128;
   Int_t           L1_SingleJet128_Prescl;
   Int_t           L1_SingleJet176;
   Int_t           L1_SingleJet176_Prescl;
   Int_t           L1_SingleJet200;
   Int_t           L1_SingleJet200_Prescl;
   Int_t           L1_SingleJet36;
   Int_t           L1_SingleJet36_Prescl;
   Int_t           L1_SingleJet52;
   Int_t           L1_SingleJet52_Prescl;
   Int_t           L1_SingleJet68;
   Int_t           L1_SingleJet68_Prescl;
   Int_t           L1_SingleJet92;
   Int_t           L1_SingleJet92_Prescl;
   Int_t           L1_SingleJetC20_NotBptxOR;
   Int_t           L1_SingleJetC20_NotBptxOR_Prescl;
   Int_t           L1_SingleJetC32_NotBptxOR;
   Int_t           L1_SingleJetC32_NotBptxOR_Prescl;
   Int_t           L1_TripleEG_14_10_8;
   Int_t           L1_TripleEG_14_10_8_Prescl;
   Int_t           L1_TripleMu0;
   Int_t           L1_TripleMu0_Prescl;
   Int_t           L1_TripleMu_5_5_3;
   Int_t           L1_TripleMu_5_5_3_Prescl;
   Int_t           L1_TripleJet_84_68_48_VBF;
   Int_t           L1_TripleJet_84_68_48_VBF_Prescl;
   Int_t           L1_TripleJet_92_76_64_VBF;
   Int_t           L1_TripleJet_92_76_64_VBF_Prescl;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;
   Int_t           L1Tech_BPTX_PreBPTX_v0;
   Int_t           L1Tech_BPTX_PreBPTX_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_v0;
   Int_t           L1Tech_BPTX_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_v0;
   Int_t           L1Tech_BPTX_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_quiet_v0;
   Int_t           L1Tech_BPTX_quiet_v0_Prescl;
   Int_t           L1Tech_CASTOR_Gap_v0;
   Int_t           L1Tech_CASTOR_Gap_v0_Prescl;
   Int_t           L1Tech_CASTOR_HaloMuon_v0;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_Prescl;
   Int_t           L1Tech_CASTOR_HighJet_v0;
   Int_t           L1Tech_CASTOR_HighJet_v0_Prescl;
   Int_t           L1Tech_CASTOR_MediumJet_v0;
   Int_t           L1Tech_CASTOR_MediumJet_v0_Prescl;
   Int_t           L1Tech_DT_GlobalOR_v0;
   Int_t           L1Tech_DT_GlobalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;
   Int_t           L1Tech_HCAL_HF_single_channel_v0;
   Int_t           L1Tech_HCAL_HF_single_channel_v0_Prescl;
   Int_t           L1Tech_HCAL_HO_totalOR_v0;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;
   Int_t           L1Tech_TOTEM_0;
   Int_t           L1Tech_TOTEM_0_Prescl;
   Int_t           L1Tech_TOTEM_1;
   Int_t           L1Tech_TOTEM_1_Prescl;
   Int_t           L1Tech_TOTEM_2;
   Int_t           L1Tech_TOTEM_2_Prescl;
   Int_t           L1Tech_TOTEM_3;
   Int_t           L1Tech_TOTEM_3_Prescl;
   Int_t           L1Tech__TTU_RB0_Cosmics_v0;
   Int_t           L1Tech__TTU_RB0_Cosmics_v0_Prescl;

   // List of branches
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_L1MuChg;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1EtTot;   //!
   TBranch        *b_L1Mht;   //!
   TBranch        *b_L1MhtPhi;   //!
   TBranch        *b_L1EtHad;   //!
   TBranch        *b_L1HfRing1EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing2EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing1EtSumNegativeEta;   //!
   TBranch        *b_L1HfRing2EtSumNegativeEta;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing1;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing1;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing2;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing2;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_digitisation_step;   //!
   TBranch        *b_digitisation_step_Prescl;   //!
   TBranch        *b_L1simulation_step;   //!
   TBranch        *b_L1simulation_step_Prescl;   //!
   TBranch        *b_digi2raw_step;   //!
   TBranch        *b_digi2raw_step_Prescl;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLTriggerFirstPath_Prescl;   //!
   TBranch        *b_HLT_Physics_v1;   //!
   TBranch        *b_HLT_Physics_v1_Prescl;   //!
   TBranch        *b_HLT_Random_v1;   //!
   TBranch        *b_HLT_Random_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_v1;   //!
   TBranch        *b_HLT_ZeroBias_v1_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_raw2digi_step;   //!
   TBranch        *b_raw2digi_step_Prescl;   //!
   TBranch        *b_L1Reco_step;   //!
   TBranch        *b_L1Reco_step_Prescl;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_AlwaysTrue_Prescl;   //!
   TBranch        *b_L1_DoubleEG6_HTT150;   //!
   TBranch        *b_L1_DoubleEG6_HTT150_Prescl;   //!
   TBranch        *b_L1_DoubleEG_15_10;   //!
   TBranch        *b_L1_DoubleEG_15_10_Prescl;   //!
   TBranch        *b_L1_DoubleEG_22_10;   //!
   TBranch        *b_L1_DoubleEG_22_10_Prescl;   //!
   TBranch        *b_L1_DoubleIsoTau28er;   //!
   TBranch        *b_L1_DoubleIsoTau28er_Prescl;   //!
   TBranch        *b_L1_DoubleIsoTau32er;   //!
   TBranch        *b_L1_DoubleIsoTau32er_Prescl;   //!
   TBranch        *b_L1_DoubleIsoTau36er;   //!
   TBranch        *b_L1_DoubleIsoTau36er_Prescl;   //!
   TBranch        *b_L1_DoubleIsoTau40er;   //!
   TBranch        *b_L1_DoubleIsoTau40er_Prescl;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_Eta1p6_WdEta18;   //!
   TBranch        *b_L1_DoubleMu0_Eta1p6_WdEta18_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_Eta1p6_WdEta18_OS;   //!
   TBranch        *b_L1_DoubleMu0_Eta1p6_WdEta18_OS_Prescl;   //!
   TBranch        *b_L1_DoubleMu6_EG6;   //!
   TBranch        *b_L1_DoubleMu6_EG6_Prescl;   //!
   TBranch        *b_L1_DoubleMu7_EG7;   //!
   TBranch        *b_L1_DoubleMu7_EG7_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen;   //!
   TBranch        *b_L1_DoubleMuOpen_Prescl;   //!
   TBranch        *b_L1_DoubleMu_10_0_WdEta18;   //!
   TBranch        *b_L1_DoubleMu_10_0_WdEta18_Prescl;   //!
   TBranch        *b_L1_DoubleMu_10_3p5;   //!
   TBranch        *b_L1_DoubleMu_10_3p5_Prescl;   //!
   TBranch        *b_L1_DoubleMu_10_Open;   //!
   TBranch        *b_L1_DoubleMu_10_Open_Prescl;   //!
   TBranch        *b_L1_DoubleMu_12_5;   //!
   TBranch        *b_L1_DoubleMu_12_5_Prescl;   //!
   TBranch        *b_L1_DoubleJetC100;   //!
   TBranch        *b_L1_DoubleJetC100_Prescl;   //!
   TBranch        *b_L1_DoubleJetC112;   //!
   TBranch        *b_L1_DoubleJetC112_Prescl;   //!
   TBranch        *b_L1_DoubleJetC52;   //!
   TBranch        *b_L1_DoubleJetC52_Prescl;   //!
   TBranch        *b_L1_DoubleJetC56_ETM60;   //!
   TBranch        *b_L1_DoubleJetC56_ETM60_Prescl;   //!
   TBranch        *b_L1_DoubleJetC84;   //!
   TBranch        *b_L1_DoubleJetC84_Prescl;   //!
   TBranch        *b_L1_DoubleTauJet40er;   //!
   TBranch        *b_L1_DoubleTauJet40er_Prescl;   //!
   TBranch        *b_L1_EG25er_HTT100;   //!
   TBranch        *b_L1_EG25er_HTT100_Prescl;   //!
   TBranch        *b_L1_ETM100;   //!
   TBranch        *b_L1_ETM100_Prescl;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM30_Prescl;   //!
   TBranch        *b_L1_ETM40;   //!
   TBranch        *b_L1_ETM40_Prescl;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM50_Prescl;   //!
   TBranch        *b_L1_ETM60;   //!
   TBranch        *b_L1_ETM60_Prescl;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM70_Prescl;   //!
   TBranch        *b_L1_ETT15_BptxAND;   //!
   TBranch        *b_L1_ETT15_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT40;   //!
   TBranch        *b_L1_ETT40_Prescl;   //!
   TBranch        *b_L1_ETT60;   //!
   TBranch        *b_L1_ETT60_Prescl;   //!
   TBranch        *b_L1_HTT100;   //!
   TBranch        *b_L1_HTT100_Prescl;   //!
   TBranch        *b_L1_HTT125;   //!
   TBranch        *b_L1_HTT125_Prescl;   //!
   TBranch        *b_L1_HTT150;   //!
   TBranch        *b_L1_HTT150_Prescl;   //!
   TBranch        *b_L1_HTT175;   //!
   TBranch        *b_L1_HTT175_Prescl;   //!
   TBranch        *b_L1_HTT200;   //!
   TBranch        *b_L1_HTT200_Prescl;   //!
   TBranch        *b_L1_HTT250;   //!
   TBranch        *b_L1_HTT250_Prescl;   //!
   TBranch        *b_L1_HTT75;   //!
   TBranch        *b_L1_HTT75_Prescl;   //!
   TBranch        *b_L1_IsoEG20er_TauJet20er_NotWdEta0;   //!
   TBranch        *b_L1_IsoEG20er_TauJet20er_NotWdEta0_Prescl;   //!
   TBranch        *b_L1_Mu0er_ETM40;   //!
   TBranch        *b_L1_Mu0er_ETM40_Prescl;   //!
   TBranch        *b_L1_Mu0er_ETM55;   //!
   TBranch        *b_L1_Mu0er_ETM55_Prescl;   //!
   TBranch        *b_L1_Mu10er_ETM30;   //!
   TBranch        *b_L1_Mu10er_ETM30_Prescl;   //!
   TBranch        *b_L1_Mu10er_ETM50;   //!
   TBranch        *b_L1_Mu10er_ETM50_Prescl;   //!
   TBranch        *b_L1_Mu12_EG10;   //!
   TBranch        *b_L1_Mu12_EG10_Prescl;   //!
   TBranch        *b_L1_Mu14er_ETM30;   //!
   TBranch        *b_L1_Mu14er_ETM30_Prescl;   //!
   TBranch        *b_L1_Mu16er_IsoTau28er;   //!
   TBranch        *b_L1_Mu16er_IsoTau28er_Prescl;   //!
   TBranch        *b_L1_Mu16er_IsoTau32er;   //!
   TBranch        *b_L1_Mu16er_IsoTau32er_Prescl;   //!
   TBranch        *b_L1_Mu16er_TauJet20er;   //!
   TBranch        *b_L1_Mu16er_TauJet20er_Prescl;   //!
   TBranch        *b_L1_Mu20_EG10;   //!
   TBranch        *b_L1_Mu20_EG10_Prescl;   //!
   TBranch        *b_L1_Mu3_JetC16_WdEtaPhi2;   //!
   TBranch        *b_L1_Mu3_JetC16_WdEtaPhi2_Prescl;   //!
   TBranch        *b_L1_Mu3_JetC52_WdEtaPhi2;   //!
   TBranch        *b_L1_Mu3_JetC52_WdEtaPhi2_Prescl;   //!
   TBranch        *b_L1_Mu4_EG18;   //!
   TBranch        *b_L1_Mu4_EG18_Prescl;   //!
   TBranch        *b_L1_Mu5_DoubleEG5;   //!
   TBranch        *b_L1_Mu5_DoubleEG5_Prescl;   //!
   TBranch        *b_L1_Mu5_EG15;   //!
   TBranch        *b_L1_Mu5_EG15_Prescl;   //!
   TBranch        *b_L1_Mu5_EG20;   //!
   TBranch        *b_L1_Mu5_EG20_Prescl;   //!
   TBranch        *b_L1_Mu5_IsoEG18;   //!
   TBranch        *b_L1_Mu5_IsoEG18_Prescl;   //!
   TBranch        *b_L1_Mu6_DoubleEG10;   //!
   TBranch        *b_L1_Mu6_DoubleEG10_Prescl;   //!
   TBranch        *b_L1_Mu6_HTT100;   //!
   TBranch        *b_L1_Mu6_HTT100_Prescl;   //!
   TBranch        *b_L1_Mu8_HTT50;   //!
   TBranch        *b_L1_Mu8_HTT50_Prescl;   //!
   TBranch        *b_L1_QuadMu0;   //!
   TBranch        *b_L1_QuadMu0_Prescl;   //!
   TBranch        *b_L1_QuadJetC36_TauJet52;   //!
   TBranch        *b_L1_QuadJetC36_TauJet52_Prescl;   //!
   TBranch        *b_L1_QuadJetC40;   //!
   TBranch        *b_L1_QuadJetC40_Prescl;   //!
   TBranch        *b_L1_QuadJetC60;   //!
   TBranch        *b_L1_QuadJetC60_Prescl;   //!
   TBranch        *b_L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1;   //!
   TBranch        *b_L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1_Prescl;   //!
   TBranch        *b_L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1;   //!
   TBranch        *b_L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1_Prescl;   //!
   TBranch        *b_L1_SingleEG10;   //!
   TBranch        *b_L1_SingleEG10_Prescl;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG15_Prescl;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG25;   //!
   TBranch        *b_L1_SingleEG25_Prescl;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG2_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG30_Prescl;   //!
   TBranch        *b_L1_SingleEG35;   //!
   TBranch        *b_L1_SingleEG35_Prescl;   //!
   TBranch        *b_L1_SingleEG40;   //!
   TBranch        *b_L1_SingleEG40_Prescl;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG5_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG18er;   //!
   TBranch        *b_L1_SingleIsoEG18er_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG20;   //!
   TBranch        *b_L1_SingleIsoEG20_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG20er;   //!
   TBranch        *b_L1_SingleIsoEG20er_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG22er;   //!
   TBranch        *b_L1_SingleIsoEG22er_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG25;   //!
   TBranch        *b_L1_SingleIsoEG25_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG25er;   //!
   TBranch        *b_L1_SingleIsoEG25er_Prescl;   //!
   TBranch        *b_L1_SingleIsoEG30er;   //!
   TBranch        *b_L1_SingleIsoEG30er_Prescl;   //!
   TBranch        *b_L1_SingleMu12;   //!
   TBranch        *b_L1_SingleMu12_Prescl;   //!
   TBranch        *b_L1_SingleMu14er;   //!
   TBranch        *b_L1_SingleMu14er_Prescl;   //!
   TBranch        *b_L1_SingleMu16;   //!
   TBranch        *b_L1_SingleMu16_Prescl;   //!
   TBranch        *b_L1_SingleMu16er;   //!
   TBranch        *b_L1_SingleMu16er_Prescl;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20_Prescl;   //!
   TBranch        *b_L1_SingleMu20er;   //!
   TBranch        *b_L1_SingleMu20er_Prescl;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleMu25_Prescl;   //!
   TBranch        *b_L1_SingleMu30;   //!
   TBranch        *b_L1_SingleMu30_Prescl;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu5_Prescl;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJet128;   //!
   TBranch        *b_L1_SingleJet128_Prescl;   //!
   TBranch        *b_L1_SingleJet176;   //!
   TBranch        *b_L1_SingleJet176_Prescl;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet200_Prescl;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet52;   //!
   TBranch        *b_L1_SingleJet52_Prescl;   //!
   TBranch        *b_L1_SingleJet68;   //!
   TBranch        *b_L1_SingleJet68_Prescl;   //!
   TBranch        *b_L1_SingleJet92;   //!
   TBranch        *b_L1_SingleJet92_Prescl;   //!
   TBranch        *b_L1_SingleJetC20_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC20_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_TripleEG_14_10_8;   //!
   TBranch        *b_L1_TripleEG_14_10_8_Prescl;   //!
   TBranch        *b_L1_TripleMu0;   //!
   TBranch        *b_L1_TripleMu0_Prescl;   //!
   TBranch        *b_L1_TripleMu_5_5_3;   //!
   TBranch        *b_L1_TripleMu_5_5_3_Prescl;   //!
   TBranch        *b_L1_TripleJet_84_68_48_VBF;   //!
   TBranch        *b_L1_TripleJet_84_68_48_VBF_Prescl;   //!
   TBranch        *b_L1_TripleJet_92_76_64_VBF;   //!
   TBranch        *b_L1_TripleJet_92_76_64_VBF_Prescl;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_Gap_v0;   //!
   TBranch        *b_L1Tech_CASTOR_Gap_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HighJet_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HighJet_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_MediumJet_v0;   //!
   TBranch        *b_L1Tech_CASTOR_MediumJet_v0_Prescl;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_0;   //!
   TBranch        *b_L1Tech_TOTEM_0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_1;   //!
   TBranch        *b_L1Tech_TOTEM_1_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_2;   //!
   TBranch        *b_L1Tech_TOTEM_2_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_3;   //!
   TBranch        *b_L1Tech_TOTEM_3_Prescl;   //!
   TBranch        *b_L1Tech__TTU_RB0_Cosmics_v0;   //!
   TBranch        *b_L1Tech__TTU_RB0_Cosmics_v0_Prescl;   //!

};


void setupHltTree(TTree *t,Hlts &tHlts,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("NL1IsolEm", &tHlts.NL1IsolEm, &tHlts.b_NL1IsolEm);
   t->SetBranchAddress("L1IsolEmEt", tHlts.L1IsolEmEt, &tHlts.b_L1IsolEmEt);
   t->SetBranchAddress("L1IsolEmE", tHlts.L1IsolEmE, &tHlts.b_L1IsolEmE);
   t->SetBranchAddress("L1IsolEmEta", tHlts.L1IsolEmEta, &tHlts.b_L1IsolEmEta);
   t->SetBranchAddress("L1IsolEmPhi", tHlts.L1IsolEmPhi, &tHlts.b_L1IsolEmPhi);
   t->SetBranchAddress("NL1NIsolEm", &tHlts.NL1NIsolEm, &tHlts.b_NL1NIsolEm);
   t->SetBranchAddress("L1NIsolEmEt", tHlts.L1NIsolEmEt, &tHlts.b_L1NIsolEmEt);
   t->SetBranchAddress("L1NIsolEmE", tHlts.L1NIsolEmE, &tHlts.b_L1NIsolEmE);
   t->SetBranchAddress("L1NIsolEmEta", tHlts.L1NIsolEmEta, &tHlts.b_L1NIsolEmEta);
   t->SetBranchAddress("L1NIsolEmPhi", tHlts.L1NIsolEmPhi, &tHlts.b_L1NIsolEmPhi);
   t->SetBranchAddress("NL1Mu", &tHlts.NL1Mu, &tHlts.b_NL1Mu);
   t->SetBranchAddress("L1MuPt", tHlts.L1MuPt, &tHlts.b_L1MuPt);
   t->SetBranchAddress("L1MuE", tHlts.L1MuE, &tHlts.b_L1MuE);
   t->SetBranchAddress("L1MuEta", tHlts.L1MuEta, &tHlts.b_L1MuEta);
   t->SetBranchAddress("L1MuPhi", tHlts.L1MuPhi, &tHlts.b_L1MuPhi);
   t->SetBranchAddress("L1MuIsol", tHlts.L1MuIsol, &tHlts.b_L1MuIsol);
   t->SetBranchAddress("L1MuMip", tHlts.L1MuMip, &tHlts.b_L1MuMip);
   t->SetBranchAddress("L1MuFor", tHlts.L1MuFor, &tHlts.b_L1MuFor);
   t->SetBranchAddress("L1MuRPC", tHlts.L1MuRPC, &tHlts.b_L1MuRPC);
   t->SetBranchAddress("L1MuQal", tHlts.L1MuQal, &tHlts.b_L1MuQal);
   t->SetBranchAddress("L1MuChg", tHlts.L1MuChg, &tHlts.b_L1MuChg);
   t->SetBranchAddress("NL1CenJet", &tHlts.NL1CenJet, &tHlts.b_NL1CenJet);
   t->SetBranchAddress("L1CenJetEt", tHlts.L1CenJetEt, &tHlts.b_L1CenJetEt);
   t->SetBranchAddress("L1CenJetE", tHlts.L1CenJetE, &tHlts.b_L1CenJetE);
   t->SetBranchAddress("L1CenJetEta", tHlts.L1CenJetEta, &tHlts.b_L1CenJetEta);
   t->SetBranchAddress("L1CenJetPhi", tHlts.L1CenJetPhi, &tHlts.b_L1CenJetPhi);
   t->SetBranchAddress("NL1ForJet", &tHlts.NL1ForJet, &tHlts.b_NL1ForJet);
   t->SetBranchAddress("L1ForJetEt", tHlts.L1ForJetEt, &tHlts.b_L1ForJetEt);
   t->SetBranchAddress("L1ForJetE", tHlts.L1ForJetE, &tHlts.b_L1ForJetE);
   t->SetBranchAddress("L1ForJetEta", tHlts.L1ForJetEta, &tHlts.b_L1ForJetEta);
   t->SetBranchAddress("L1ForJetPhi", tHlts.L1ForJetPhi, &tHlts.b_L1ForJetPhi);
   t->SetBranchAddress("NL1Tau", &tHlts.NL1Tau, &tHlts.b_NL1Tau);
   t->SetBranchAddress("L1TauEt", tHlts.L1TauEt, &tHlts.b_L1TauEt);
   t->SetBranchAddress("L1TauE", tHlts.L1TauE, &tHlts.b_L1TauE);
   t->SetBranchAddress("L1TauEta", tHlts.L1TauEta, &tHlts.b_L1TauEta);
   t->SetBranchAddress("L1TauPhi", tHlts.L1TauPhi, &tHlts.b_L1TauPhi);
   t->SetBranchAddress("L1Met", &tHlts.L1Met, &tHlts.b_L1Met);
   t->SetBranchAddress("L1MetPhi", &tHlts.L1MetPhi, &tHlts.b_L1MetPhi);
   t->SetBranchAddress("L1EtTot", &tHlts.L1EtTot, &tHlts.b_L1EtTot);
   t->SetBranchAddress("L1Mht", &tHlts.L1Mht, &tHlts.b_L1Mht);
   t->SetBranchAddress("L1MhtPhi", &tHlts.L1MhtPhi, &tHlts.b_L1MhtPhi);
   t->SetBranchAddress("L1EtHad", &tHlts.L1EtHad, &tHlts.b_L1EtHad);
   t->SetBranchAddress("L1HfRing1EtSumPositiveEta", &tHlts.L1HfRing1EtSumPositiveEta, &tHlts.b_L1HfRing1EtSumPositiveEta);
   t->SetBranchAddress("L1HfRing2EtSumPositiveEta", &tHlts.L1HfRing2EtSumPositiveEta, &tHlts.b_L1HfRing2EtSumPositiveEta);
   t->SetBranchAddress("L1HfRing1EtSumNegativeEta", &tHlts.L1HfRing1EtSumNegativeEta, &tHlts.b_L1HfRing1EtSumNegativeEta);
   t->SetBranchAddress("L1HfRing2EtSumNegativeEta", &tHlts.L1HfRing2EtSumNegativeEta, &tHlts.b_L1HfRing2EtSumNegativeEta);
   t->SetBranchAddress("L1HfTowerCountPositiveEtaRing1", &tHlts.L1HfTowerCountPositiveEtaRing1, &tHlts.b_L1HfTowerCountPositiveEtaRing1);
   t->SetBranchAddress("L1HfTowerCountNegativeEtaRing1", &tHlts.L1HfTowerCountNegativeEtaRing1, &tHlts.b_L1HfTowerCountNegativeEtaRing1);
   t->SetBranchAddress("L1HfTowerCountPositiveEtaRing2", &tHlts.L1HfTowerCountPositiveEtaRing2, &tHlts.b_L1HfTowerCountPositiveEtaRing2);
   t->SetBranchAddress("L1HfTowerCountNegativeEtaRing2", &tHlts.L1HfTowerCountNegativeEtaRing2, &tHlts.b_L1HfTowerCountNegativeEtaRing2);
   t->SetBranchAddress("Event", &tHlts.Event, &tHlts.b_Event);
   t->SetBranchAddress("LumiBlock", &tHlts.LumiBlock, &tHlts.b_LumiBlock);
   t->SetBranchAddress("Run", &tHlts.Run, &tHlts.b_Run);
   t->SetBranchAddress("Bx", &tHlts.Bx, &tHlts.b_Bx);
   t->SetBranchAddress("Orbit", &tHlts.Orbit, &tHlts.b_Orbit);
   t->SetBranchAddress("AvgInstDelLumi", &tHlts.AvgInstDelLumi, &tHlts.b_AvgInstDelLumi);
   t->SetBranchAddress("digitisation_step", &tHlts.digitisation_step, &tHlts.b_digitisation_step);
   t->SetBranchAddress("digitisation_step_Prescl", &tHlts.digitisation_step_Prescl, &tHlts.b_digitisation_step_Prescl);
   t->SetBranchAddress("L1simulation_step", &tHlts.L1simulation_step, &tHlts.b_L1simulation_step);
   t->SetBranchAddress("L1simulation_step_Prescl", &tHlts.L1simulation_step_Prescl, &tHlts.b_L1simulation_step_Prescl);
   t->SetBranchAddress("digi2raw_step", &tHlts.digi2raw_step, &tHlts.b_digi2raw_step);
   t->SetBranchAddress("digi2raw_step_Prescl", &tHlts.digi2raw_step_Prescl, &tHlts.b_digi2raw_step_Prescl);
   t->SetBranchAddress("HLTriggerFirstPath", &tHlts.HLTriggerFirstPath, &tHlts.b_HLTriggerFirstPath);
   t->SetBranchAddress("HLTriggerFirstPath_Prescl", &tHlts.HLTriggerFirstPath_Prescl, &tHlts.b_HLTriggerFirstPath_Prescl);
   t->SetBranchAddress("HLT_Physics_v1", &tHlts.HLT_Physics_v1, &tHlts.b_HLT_Physics_v1);
   t->SetBranchAddress("HLT_Physics_v1_Prescl", &tHlts.HLT_Physics_v1_Prescl, &tHlts.b_HLT_Physics_v1_Prescl);
   t->SetBranchAddress("HLT_Random_v1", &tHlts.HLT_Random_v1, &tHlts.b_HLT_Random_v1);
   t->SetBranchAddress("HLT_Random_v1_Prescl", &tHlts.HLT_Random_v1_Prescl, &tHlts.b_HLT_Random_v1_Prescl);
   t->SetBranchAddress("HLT_ZeroBias_v1", &tHlts.HLT_ZeroBias_v1, &tHlts.b_HLT_ZeroBias_v1);
   t->SetBranchAddress("HLT_ZeroBias_v1_Prescl", &tHlts.HLT_ZeroBias_v1_Prescl, &tHlts.b_HLT_ZeroBias_v1_Prescl);
   t->SetBranchAddress("HLTriggerFinalPath", &tHlts.HLTriggerFinalPath, &tHlts.b_HLTriggerFinalPath);
   t->SetBranchAddress("HLTriggerFinalPath_Prescl", &tHlts.HLTriggerFinalPath_Prescl, &tHlts.b_HLTriggerFinalPath_Prescl);
   t->SetBranchAddress("raw2digi_step", &tHlts.raw2digi_step, &tHlts.b_raw2digi_step);
   t->SetBranchAddress("raw2digi_step_Prescl", &tHlts.raw2digi_step_Prescl, &tHlts.b_raw2digi_step_Prescl);
   t->SetBranchAddress("L1Reco_step", &tHlts.L1Reco_step, &tHlts.b_L1Reco_step);
   t->SetBranchAddress("L1Reco_step_Prescl", &tHlts.L1Reco_step_Prescl, &tHlts.b_L1Reco_step_Prescl);
   t->SetBranchAddress("L1_AlwaysTrue", &tHlts.L1_AlwaysTrue, &tHlts.b_L1_AlwaysTrue);
   t->SetBranchAddress("L1_AlwaysTrue_Prescl", &tHlts.L1_AlwaysTrue_Prescl, &tHlts.b_L1_AlwaysTrue_Prescl);
   t->SetBranchAddress("L1_DoubleEG6_HTT150", &tHlts.L1_DoubleEG6_HTT150, &tHlts.b_L1_DoubleEG6_HTT150);
   t->SetBranchAddress("L1_DoubleEG6_HTT150_Prescl", &tHlts.L1_DoubleEG6_HTT150_Prescl, &tHlts.b_L1_DoubleEG6_HTT150_Prescl);
   t->SetBranchAddress("L1_DoubleEG_15_10", &tHlts.L1_DoubleEG_15_10, &tHlts.b_L1_DoubleEG_15_10);
   t->SetBranchAddress("L1_DoubleEG_15_10_Prescl", &tHlts.L1_DoubleEG_15_10_Prescl, &tHlts.b_L1_DoubleEG_15_10_Prescl);
   t->SetBranchAddress("L1_DoubleEG_22_10", &tHlts.L1_DoubleEG_22_10, &tHlts.b_L1_DoubleEG_22_10);
   t->SetBranchAddress("L1_DoubleEG_22_10_Prescl", &tHlts.L1_DoubleEG_22_10_Prescl, &tHlts.b_L1_DoubleEG_22_10_Prescl);
   t->SetBranchAddress("L1_DoubleIsoTau28er", &tHlts.L1_DoubleIsoTau28er, &tHlts.b_L1_DoubleIsoTau28er);
   t->SetBranchAddress("L1_DoubleIsoTau28er_Prescl", &tHlts.L1_DoubleIsoTau28er_Prescl, &tHlts.b_L1_DoubleIsoTau28er_Prescl);
   t->SetBranchAddress("L1_DoubleIsoTau32er", &tHlts.L1_DoubleIsoTau32er, &tHlts.b_L1_DoubleIsoTau32er);
   t->SetBranchAddress("L1_DoubleIsoTau32er_Prescl", &tHlts.L1_DoubleIsoTau32er_Prescl, &tHlts.b_L1_DoubleIsoTau32er_Prescl);
   t->SetBranchAddress("L1_DoubleIsoTau36er", &tHlts.L1_DoubleIsoTau36er, &tHlts.b_L1_DoubleIsoTau36er);
   t->SetBranchAddress("L1_DoubleIsoTau36er_Prescl", &tHlts.L1_DoubleIsoTau36er_Prescl, &tHlts.b_L1_DoubleIsoTau36er_Prescl);
   t->SetBranchAddress("L1_DoubleIsoTau40er", &tHlts.L1_DoubleIsoTau40er, &tHlts.b_L1_DoubleIsoTau40er);
   t->SetBranchAddress("L1_DoubleIsoTau40er_Prescl", &tHlts.L1_DoubleIsoTau40er_Prescl, &tHlts.b_L1_DoubleIsoTau40er_Prescl);
   t->SetBranchAddress("L1_DoubleMu0", &tHlts.L1_DoubleMu0, &tHlts.b_L1_DoubleMu0);
   t->SetBranchAddress("L1_DoubleMu0_Prescl", &tHlts.L1_DoubleMu0_Prescl, &tHlts.b_L1_DoubleMu0_Prescl);
   t->SetBranchAddress("L1_DoubleMu0_Eta1p6_WdEta18", &tHlts.L1_DoubleMu0_Eta1p6_WdEta18, &tHlts.b_L1_DoubleMu0_Eta1p6_WdEta18);
   t->SetBranchAddress("L1_DoubleMu0_Eta1p6_WdEta18_Prescl", &tHlts.L1_DoubleMu0_Eta1p6_WdEta18_Prescl, &tHlts.b_L1_DoubleMu0_Eta1p6_WdEta18_Prescl);
   t->SetBranchAddress("L1_DoubleMu0_Eta1p6_WdEta18_OS", &tHlts.L1_DoubleMu0_Eta1p6_WdEta18_OS, &tHlts.b_L1_DoubleMu0_Eta1p6_WdEta18_OS);
   t->SetBranchAddress("L1_DoubleMu0_Eta1p6_WdEta18_OS_Prescl", &tHlts.L1_DoubleMu0_Eta1p6_WdEta18_OS_Prescl, &tHlts.b_L1_DoubleMu0_Eta1p6_WdEta18_OS_Prescl);
   t->SetBranchAddress("L1_DoubleMu6_EG6", &tHlts.L1_DoubleMu6_EG6, &tHlts.b_L1_DoubleMu6_EG6);
   t->SetBranchAddress("L1_DoubleMu6_EG6_Prescl", &tHlts.L1_DoubleMu6_EG6_Prescl, &tHlts.b_L1_DoubleMu6_EG6_Prescl);
   t->SetBranchAddress("L1_DoubleMu7_EG7", &tHlts.L1_DoubleMu7_EG7, &tHlts.b_L1_DoubleMu7_EG7);
   t->SetBranchAddress("L1_DoubleMu7_EG7_Prescl", &tHlts.L1_DoubleMu7_EG7_Prescl, &tHlts.b_L1_DoubleMu7_EG7_Prescl);
   t->SetBranchAddress("L1_DoubleMuOpen", &tHlts.L1_DoubleMuOpen, &tHlts.b_L1_DoubleMuOpen);
   t->SetBranchAddress("L1_DoubleMuOpen_Prescl", &tHlts.L1_DoubleMuOpen_Prescl, &tHlts.b_L1_DoubleMuOpen_Prescl);
   t->SetBranchAddress("L1_DoubleMu_10_0_WdEta18", &tHlts.L1_DoubleMu_10_0_WdEta18, &tHlts.b_L1_DoubleMu_10_0_WdEta18);
   t->SetBranchAddress("L1_DoubleMu_10_0_WdEta18_Prescl", &tHlts.L1_DoubleMu_10_0_WdEta18_Prescl, &tHlts.b_L1_DoubleMu_10_0_WdEta18_Prescl);
   t->SetBranchAddress("L1_DoubleMu_10_3p5", &tHlts.L1_DoubleMu_10_3p5, &tHlts.b_L1_DoubleMu_10_3p5);
   t->SetBranchAddress("L1_DoubleMu_10_3p5_Prescl", &tHlts.L1_DoubleMu_10_3p5_Prescl, &tHlts.b_L1_DoubleMu_10_3p5_Prescl);
   t->SetBranchAddress("L1_DoubleMu_10_Open", &tHlts.L1_DoubleMu_10_Open, &tHlts.b_L1_DoubleMu_10_Open);
   t->SetBranchAddress("L1_DoubleMu_10_Open_Prescl", &tHlts.L1_DoubleMu_10_Open_Prescl, &tHlts.b_L1_DoubleMu_10_Open_Prescl);
   t->SetBranchAddress("L1_DoubleMu_12_5", &tHlts.L1_DoubleMu_12_5, &tHlts.b_L1_DoubleMu_12_5);
   t->SetBranchAddress("L1_DoubleMu_12_5_Prescl", &tHlts.L1_DoubleMu_12_5_Prescl, &tHlts.b_L1_DoubleMu_12_5_Prescl);
   t->SetBranchAddress("L1_DoubleJetC100", &tHlts.L1_DoubleJetC100, &tHlts.b_L1_DoubleJetC100);
   t->SetBranchAddress("L1_DoubleJetC100_Prescl", &tHlts.L1_DoubleJetC100_Prescl, &tHlts.b_L1_DoubleJetC100_Prescl);
   t->SetBranchAddress("L1_DoubleJetC112", &tHlts.L1_DoubleJetC112, &tHlts.b_L1_DoubleJetC112);
   t->SetBranchAddress("L1_DoubleJetC112_Prescl", &tHlts.L1_DoubleJetC112_Prescl, &tHlts.b_L1_DoubleJetC112_Prescl);
   t->SetBranchAddress("L1_DoubleJetC52", &tHlts.L1_DoubleJetC52, &tHlts.b_L1_DoubleJetC52);
   t->SetBranchAddress("L1_DoubleJetC52_Prescl", &tHlts.L1_DoubleJetC52_Prescl, &tHlts.b_L1_DoubleJetC52_Prescl);
   t->SetBranchAddress("L1_DoubleJetC56_ETM60", &tHlts.L1_DoubleJetC56_ETM60, &tHlts.b_L1_DoubleJetC56_ETM60);
   t->SetBranchAddress("L1_DoubleJetC56_ETM60_Prescl", &tHlts.L1_DoubleJetC56_ETM60_Prescl, &tHlts.b_L1_DoubleJetC56_ETM60_Prescl);
   t->SetBranchAddress("L1_DoubleJetC84", &tHlts.L1_DoubleJetC84, &tHlts.b_L1_DoubleJetC84);
   t->SetBranchAddress("L1_DoubleJetC84_Prescl", &tHlts.L1_DoubleJetC84_Prescl, &tHlts.b_L1_DoubleJetC84_Prescl);
   t->SetBranchAddress("L1_DoubleTauJet40er", &tHlts.L1_DoubleTauJet40er, &tHlts.b_L1_DoubleTauJet40er);
   t->SetBranchAddress("L1_DoubleTauJet40er_Prescl", &tHlts.L1_DoubleTauJet40er_Prescl, &tHlts.b_L1_DoubleTauJet40er_Prescl);
   t->SetBranchAddress("L1_EG25er_HTT100", &tHlts.L1_EG25er_HTT100, &tHlts.b_L1_EG25er_HTT100);
   t->SetBranchAddress("L1_EG25er_HTT100_Prescl", &tHlts.L1_EG25er_HTT100_Prescl, &tHlts.b_L1_EG25er_HTT100_Prescl);
   t->SetBranchAddress("L1_ETM100", &tHlts.L1_ETM100, &tHlts.b_L1_ETM100);
   t->SetBranchAddress("L1_ETM100_Prescl", &tHlts.L1_ETM100_Prescl, &tHlts.b_L1_ETM100_Prescl);
   t->SetBranchAddress("L1_ETM30", &tHlts.L1_ETM30, &tHlts.b_L1_ETM30);
   t->SetBranchAddress("L1_ETM30_Prescl", &tHlts.L1_ETM30_Prescl, &tHlts.b_L1_ETM30_Prescl);
   t->SetBranchAddress("L1_ETM40", &tHlts.L1_ETM40, &tHlts.b_L1_ETM40);
   t->SetBranchAddress("L1_ETM40_Prescl", &tHlts.L1_ETM40_Prescl, &tHlts.b_L1_ETM40_Prescl);
   t->SetBranchAddress("L1_ETM50", &tHlts.L1_ETM50, &tHlts.b_L1_ETM50);
   t->SetBranchAddress("L1_ETM50_Prescl", &tHlts.L1_ETM50_Prescl, &tHlts.b_L1_ETM50_Prescl);
   t->SetBranchAddress("L1_ETM60", &tHlts.L1_ETM60, &tHlts.b_L1_ETM60);
   t->SetBranchAddress("L1_ETM60_Prescl", &tHlts.L1_ETM60_Prescl, &tHlts.b_L1_ETM60_Prescl);
   t->SetBranchAddress("L1_ETM70", &tHlts.L1_ETM70, &tHlts.b_L1_ETM70);
   t->SetBranchAddress("L1_ETM70_Prescl", &tHlts.L1_ETM70_Prescl, &tHlts.b_L1_ETM70_Prescl);
   t->SetBranchAddress("L1_ETT15_BptxAND", &tHlts.L1_ETT15_BptxAND, &tHlts.b_L1_ETT15_BptxAND);
   t->SetBranchAddress("L1_ETT15_BptxAND_Prescl", &tHlts.L1_ETT15_BptxAND_Prescl, &tHlts.b_L1_ETT15_BptxAND_Prescl);
   t->SetBranchAddress("L1_ETT40", &tHlts.L1_ETT40, &tHlts.b_L1_ETT40);
   t->SetBranchAddress("L1_ETT40_Prescl", &tHlts.L1_ETT40_Prescl, &tHlts.b_L1_ETT40_Prescl);
   t->SetBranchAddress("L1_ETT60", &tHlts.L1_ETT60, &tHlts.b_L1_ETT60);
   t->SetBranchAddress("L1_ETT60_Prescl", &tHlts.L1_ETT60_Prescl, &tHlts.b_L1_ETT60_Prescl);
   t->SetBranchAddress("L1_HTT100", &tHlts.L1_HTT100, &tHlts.b_L1_HTT100);
   t->SetBranchAddress("L1_HTT100_Prescl", &tHlts.L1_HTT100_Prescl, &tHlts.b_L1_HTT100_Prescl);
   t->SetBranchAddress("L1_HTT125", &tHlts.L1_HTT125, &tHlts.b_L1_HTT125);
   t->SetBranchAddress("L1_HTT125_Prescl", &tHlts.L1_HTT125_Prescl, &tHlts.b_L1_HTT125_Prescl);
   t->SetBranchAddress("L1_HTT150", &tHlts.L1_HTT150, &tHlts.b_L1_HTT150);
   t->SetBranchAddress("L1_HTT150_Prescl", &tHlts.L1_HTT150_Prescl, &tHlts.b_L1_HTT150_Prescl);
   t->SetBranchAddress("L1_HTT175", &tHlts.L1_HTT175, &tHlts.b_L1_HTT175);
   t->SetBranchAddress("L1_HTT175_Prescl", &tHlts.L1_HTT175_Prescl, &tHlts.b_L1_HTT175_Prescl);
   t->SetBranchAddress("L1_HTT200", &tHlts.L1_HTT200, &tHlts.b_L1_HTT200);
   t->SetBranchAddress("L1_HTT200_Prescl", &tHlts.L1_HTT200_Prescl, &tHlts.b_L1_HTT200_Prescl);
   t->SetBranchAddress("L1_HTT250", &tHlts.L1_HTT250, &tHlts.b_L1_HTT250);
   t->SetBranchAddress("L1_HTT250_Prescl", &tHlts.L1_HTT250_Prescl, &tHlts.b_L1_HTT250_Prescl);
   t->SetBranchAddress("L1_HTT75", &tHlts.L1_HTT75, &tHlts.b_L1_HTT75);
   t->SetBranchAddress("L1_HTT75_Prescl", &tHlts.L1_HTT75_Prescl, &tHlts.b_L1_HTT75_Prescl);
   t->SetBranchAddress("L1_IsoEG20er_TauJet20er_NotWdEta0", &tHlts.L1_IsoEG20er_TauJet20er_NotWdEta0, &tHlts.b_L1_IsoEG20er_TauJet20er_NotWdEta0);
   t->SetBranchAddress("L1_IsoEG20er_TauJet20er_NotWdEta0_Prescl", &tHlts.L1_IsoEG20er_TauJet20er_NotWdEta0_Prescl, &tHlts.b_L1_IsoEG20er_TauJet20er_NotWdEta0_Prescl);
   t->SetBranchAddress("L1_Mu0er_ETM40", &tHlts.L1_Mu0er_ETM40, &tHlts.b_L1_Mu0er_ETM40);
   t->SetBranchAddress("L1_Mu0er_ETM40_Prescl", &tHlts.L1_Mu0er_ETM40_Prescl, &tHlts.b_L1_Mu0er_ETM40_Prescl);
   t->SetBranchAddress("L1_Mu0er_ETM55", &tHlts.L1_Mu0er_ETM55, &tHlts.b_L1_Mu0er_ETM55);
   t->SetBranchAddress("L1_Mu0er_ETM55_Prescl", &tHlts.L1_Mu0er_ETM55_Prescl, &tHlts.b_L1_Mu0er_ETM55_Prescl);
   t->SetBranchAddress("L1_Mu10er_ETM30", &tHlts.L1_Mu10er_ETM30, &tHlts.b_L1_Mu10er_ETM30);
   t->SetBranchAddress("L1_Mu10er_ETM30_Prescl", &tHlts.L1_Mu10er_ETM30_Prescl, &tHlts.b_L1_Mu10er_ETM30_Prescl);
   t->SetBranchAddress("L1_Mu10er_ETM50", &tHlts.L1_Mu10er_ETM50, &tHlts.b_L1_Mu10er_ETM50);
   t->SetBranchAddress("L1_Mu10er_ETM50_Prescl", &tHlts.L1_Mu10er_ETM50_Prescl, &tHlts.b_L1_Mu10er_ETM50_Prescl);
   t->SetBranchAddress("L1_Mu12_EG10", &tHlts.L1_Mu12_EG10, &tHlts.b_L1_Mu12_EG10);
   t->SetBranchAddress("L1_Mu12_EG10_Prescl", &tHlts.L1_Mu12_EG10_Prescl, &tHlts.b_L1_Mu12_EG10_Prescl);
   t->SetBranchAddress("L1_Mu14er_ETM30", &tHlts.L1_Mu14er_ETM30, &tHlts.b_L1_Mu14er_ETM30);
   t->SetBranchAddress("L1_Mu14er_ETM30_Prescl", &tHlts.L1_Mu14er_ETM30_Prescl, &tHlts.b_L1_Mu14er_ETM30_Prescl);
   t->SetBranchAddress("L1_Mu16er_IsoTau28er", &tHlts.L1_Mu16er_IsoTau28er, &tHlts.b_L1_Mu16er_IsoTau28er);
   t->SetBranchAddress("L1_Mu16er_IsoTau28er_Prescl", &tHlts.L1_Mu16er_IsoTau28er_Prescl, &tHlts.b_L1_Mu16er_IsoTau28er_Prescl);
   t->SetBranchAddress("L1_Mu16er_IsoTau32er", &tHlts.L1_Mu16er_IsoTau32er, &tHlts.b_L1_Mu16er_IsoTau32er);
   t->SetBranchAddress("L1_Mu16er_IsoTau32er_Prescl", &tHlts.L1_Mu16er_IsoTau32er_Prescl, &tHlts.b_L1_Mu16er_IsoTau32er_Prescl);
   t->SetBranchAddress("L1_Mu16er_TauJet20er", &tHlts.L1_Mu16er_TauJet20er, &tHlts.b_L1_Mu16er_TauJet20er);
   t->SetBranchAddress("L1_Mu16er_TauJet20er_Prescl", &tHlts.L1_Mu16er_TauJet20er_Prescl, &tHlts.b_L1_Mu16er_TauJet20er_Prescl);
   t->SetBranchAddress("L1_Mu20_EG10", &tHlts.L1_Mu20_EG10, &tHlts.b_L1_Mu20_EG10);
   t->SetBranchAddress("L1_Mu20_EG10_Prescl", &tHlts.L1_Mu20_EG10_Prescl, &tHlts.b_L1_Mu20_EG10_Prescl);
   t->SetBranchAddress("L1_Mu3_JetC16_WdEtaPhi2", &tHlts.L1_Mu3_JetC16_WdEtaPhi2, &tHlts.b_L1_Mu3_JetC16_WdEtaPhi2);
   t->SetBranchAddress("L1_Mu3_JetC16_WdEtaPhi2_Prescl", &tHlts.L1_Mu3_JetC16_WdEtaPhi2_Prescl, &tHlts.b_L1_Mu3_JetC16_WdEtaPhi2_Prescl);
   t->SetBranchAddress("L1_Mu3_JetC52_WdEtaPhi2", &tHlts.L1_Mu3_JetC52_WdEtaPhi2, &tHlts.b_L1_Mu3_JetC52_WdEtaPhi2);
   t->SetBranchAddress("L1_Mu3_JetC52_WdEtaPhi2_Prescl", &tHlts.L1_Mu3_JetC52_WdEtaPhi2_Prescl, &tHlts.b_L1_Mu3_JetC52_WdEtaPhi2_Prescl);
   t->SetBranchAddress("L1_Mu4_EG18", &tHlts.L1_Mu4_EG18, &tHlts.b_L1_Mu4_EG18);
   t->SetBranchAddress("L1_Mu4_EG18_Prescl", &tHlts.L1_Mu4_EG18_Prescl, &tHlts.b_L1_Mu4_EG18_Prescl);
   t->SetBranchAddress("L1_Mu5_DoubleEG5", &tHlts.L1_Mu5_DoubleEG5, &tHlts.b_L1_Mu5_DoubleEG5);
   t->SetBranchAddress("L1_Mu5_DoubleEG5_Prescl", &tHlts.L1_Mu5_DoubleEG5_Prescl, &tHlts.b_L1_Mu5_DoubleEG5_Prescl);
   t->SetBranchAddress("L1_Mu5_EG15", &tHlts.L1_Mu5_EG15, &tHlts.b_L1_Mu5_EG15);
   t->SetBranchAddress("L1_Mu5_EG15_Prescl", &tHlts.L1_Mu5_EG15_Prescl, &tHlts.b_L1_Mu5_EG15_Prescl);
   t->SetBranchAddress("L1_Mu5_EG20", &tHlts.L1_Mu5_EG20, &tHlts.b_L1_Mu5_EG20);
   t->SetBranchAddress("L1_Mu5_EG20_Prescl", &tHlts.L1_Mu5_EG20_Prescl, &tHlts.b_L1_Mu5_EG20_Prescl);
   t->SetBranchAddress("L1_Mu5_IsoEG18", &tHlts.L1_Mu5_IsoEG18, &tHlts.b_L1_Mu5_IsoEG18);
   t->SetBranchAddress("L1_Mu5_IsoEG18_Prescl", &tHlts.L1_Mu5_IsoEG18_Prescl, &tHlts.b_L1_Mu5_IsoEG18_Prescl);
   t->SetBranchAddress("L1_Mu6_DoubleEG10", &tHlts.L1_Mu6_DoubleEG10, &tHlts.b_L1_Mu6_DoubleEG10);
   t->SetBranchAddress("L1_Mu6_DoubleEG10_Prescl", &tHlts.L1_Mu6_DoubleEG10_Prescl, &tHlts.b_L1_Mu6_DoubleEG10_Prescl);
   t->SetBranchAddress("L1_Mu6_HTT100", &tHlts.L1_Mu6_HTT100, &tHlts.b_L1_Mu6_HTT100);
   t->SetBranchAddress("L1_Mu6_HTT100_Prescl", &tHlts.L1_Mu6_HTT100_Prescl, &tHlts.b_L1_Mu6_HTT100_Prescl);
   t->SetBranchAddress("L1_Mu8_HTT50", &tHlts.L1_Mu8_HTT50, &tHlts.b_L1_Mu8_HTT50);
   t->SetBranchAddress("L1_Mu8_HTT50_Prescl", &tHlts.L1_Mu8_HTT50_Prescl, &tHlts.b_L1_Mu8_HTT50_Prescl);
   t->SetBranchAddress("L1_QuadMu0", &tHlts.L1_QuadMu0, &tHlts.b_L1_QuadMu0);
   t->SetBranchAddress("L1_QuadMu0_Prescl", &tHlts.L1_QuadMu0_Prescl, &tHlts.b_L1_QuadMu0_Prescl);
   t->SetBranchAddress("L1_QuadJetC36_TauJet52", &tHlts.L1_QuadJetC36_TauJet52, &tHlts.b_L1_QuadJetC36_TauJet52);
   t->SetBranchAddress("L1_QuadJetC36_TauJet52_Prescl", &tHlts.L1_QuadJetC36_TauJet52_Prescl, &tHlts.b_L1_QuadJetC36_TauJet52_Prescl);
   t->SetBranchAddress("L1_QuadJetC40", &tHlts.L1_QuadJetC40, &tHlts.b_L1_QuadJetC40);
   t->SetBranchAddress("L1_QuadJetC40_Prescl", &tHlts.L1_QuadJetC40_Prescl, &tHlts.b_L1_QuadJetC40_Prescl);
   t->SetBranchAddress("L1_QuadJetC60", &tHlts.L1_QuadJetC60, &tHlts.b_L1_QuadJetC60);
   t->SetBranchAddress("L1_QuadJetC60_Prescl", &tHlts.L1_QuadJetC60_Prescl, &tHlts.b_L1_QuadJetC60_Prescl);
   t->SetBranchAddress("L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1", &tHlts.L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1, &tHlts.b_L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1);
   t->SetBranchAddress("L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1_Prescl", &tHlts.L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1_Prescl, &tHlts.b_L1_Jet32_DoubleMu_Open_10_MuMuNotWdPhi23_JetMuWdPhi1_Prescl);
   t->SetBranchAddress("L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1", &tHlts.L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1, &tHlts.b_L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1);
   t->SetBranchAddress("L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1_Prescl", &tHlts.L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1_Prescl, &tHlts.b_L1_Jet32_MuOpen_EG10_MuEGNotWdPhi3_JetMuWdPhi1_Prescl);
   t->SetBranchAddress("L1_SingleEG10", &tHlts.L1_SingleEG10, &tHlts.b_L1_SingleEG10);
   t->SetBranchAddress("L1_SingleEG10_Prescl", &tHlts.L1_SingleEG10_Prescl, &tHlts.b_L1_SingleEG10_Prescl);
   t->SetBranchAddress("L1_SingleEG15", &tHlts.L1_SingleEG15, &tHlts.b_L1_SingleEG15);
   t->SetBranchAddress("L1_SingleEG15_Prescl", &tHlts.L1_SingleEG15_Prescl, &tHlts.b_L1_SingleEG15_Prescl);
   t->SetBranchAddress("L1_SingleEG20", &tHlts.L1_SingleEG20, &tHlts.b_L1_SingleEG20);
   t->SetBranchAddress("L1_SingleEG20_Prescl", &tHlts.L1_SingleEG20_Prescl, &tHlts.b_L1_SingleEG20_Prescl);
   t->SetBranchAddress("L1_SingleEG25", &tHlts.L1_SingleEG25, &tHlts.b_L1_SingleEG25);
   t->SetBranchAddress("L1_SingleEG25_Prescl", &tHlts.L1_SingleEG25_Prescl, &tHlts.b_L1_SingleEG25_Prescl);
   t->SetBranchAddress("L1_SingleEG2_BptxAND", &tHlts.L1_SingleEG2_BptxAND, &tHlts.b_L1_SingleEG2_BptxAND);
   t->SetBranchAddress("L1_SingleEG2_BptxAND_Prescl", &tHlts.L1_SingleEG2_BptxAND_Prescl, &tHlts.b_L1_SingleEG2_BptxAND_Prescl);
   t->SetBranchAddress("L1_SingleEG30", &tHlts.L1_SingleEG30, &tHlts.b_L1_SingleEG30);
   t->SetBranchAddress("L1_SingleEG30_Prescl", &tHlts.L1_SingleEG30_Prescl, &tHlts.b_L1_SingleEG30_Prescl);
   t->SetBranchAddress("L1_SingleEG35", &tHlts.L1_SingleEG35, &tHlts.b_L1_SingleEG35);
   t->SetBranchAddress("L1_SingleEG35_Prescl", &tHlts.L1_SingleEG35_Prescl, &tHlts.b_L1_SingleEG35_Prescl);
   t->SetBranchAddress("L1_SingleEG40", &tHlts.L1_SingleEG40, &tHlts.b_L1_SingleEG40);
   t->SetBranchAddress("L1_SingleEG40_Prescl", &tHlts.L1_SingleEG40_Prescl, &tHlts.b_L1_SingleEG40_Prescl);
   t->SetBranchAddress("L1_SingleEG5", &tHlts.L1_SingleEG5, &tHlts.b_L1_SingleEG5);
   t->SetBranchAddress("L1_SingleEG5_Prescl", &tHlts.L1_SingleEG5_Prescl, &tHlts.b_L1_SingleEG5_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG18er", &tHlts.L1_SingleIsoEG18er, &tHlts.b_L1_SingleIsoEG18er);
   t->SetBranchAddress("L1_SingleIsoEG18er_Prescl", &tHlts.L1_SingleIsoEG18er_Prescl, &tHlts.b_L1_SingleIsoEG18er_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG20", &tHlts.L1_SingleIsoEG20, &tHlts.b_L1_SingleIsoEG20);
   t->SetBranchAddress("L1_SingleIsoEG20_Prescl", &tHlts.L1_SingleIsoEG20_Prescl, &tHlts.b_L1_SingleIsoEG20_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG20er", &tHlts.L1_SingleIsoEG20er, &tHlts.b_L1_SingleIsoEG20er);
   t->SetBranchAddress("L1_SingleIsoEG20er_Prescl", &tHlts.L1_SingleIsoEG20er_Prescl, &tHlts.b_L1_SingleIsoEG20er_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG22er", &tHlts.L1_SingleIsoEG22er, &tHlts.b_L1_SingleIsoEG22er);
   t->SetBranchAddress("L1_SingleIsoEG22er_Prescl", &tHlts.L1_SingleIsoEG22er_Prescl, &tHlts.b_L1_SingleIsoEG22er_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG25", &tHlts.L1_SingleIsoEG25, &tHlts.b_L1_SingleIsoEG25);
   t->SetBranchAddress("L1_SingleIsoEG25_Prescl", &tHlts.L1_SingleIsoEG25_Prescl, &tHlts.b_L1_SingleIsoEG25_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG25er", &tHlts.L1_SingleIsoEG25er, &tHlts.b_L1_SingleIsoEG25er);
   t->SetBranchAddress("L1_SingleIsoEG25er_Prescl", &tHlts.L1_SingleIsoEG25er_Prescl, &tHlts.b_L1_SingleIsoEG25er_Prescl);
   t->SetBranchAddress("L1_SingleIsoEG30er", &tHlts.L1_SingleIsoEG30er, &tHlts.b_L1_SingleIsoEG30er);
   t->SetBranchAddress("L1_SingleIsoEG30er_Prescl", &tHlts.L1_SingleIsoEG30er_Prescl, &tHlts.b_L1_SingleIsoEG30er_Prescl);
   t->SetBranchAddress("L1_SingleMu12", &tHlts.L1_SingleMu12, &tHlts.b_L1_SingleMu12);
   t->SetBranchAddress("L1_SingleMu12_Prescl", &tHlts.L1_SingleMu12_Prescl, &tHlts.b_L1_SingleMu12_Prescl);
   t->SetBranchAddress("L1_SingleMu14er", &tHlts.L1_SingleMu14er, &tHlts.b_L1_SingleMu14er);
   t->SetBranchAddress("L1_SingleMu14er_Prescl", &tHlts.L1_SingleMu14er_Prescl, &tHlts.b_L1_SingleMu14er_Prescl);
   t->SetBranchAddress("L1_SingleMu16", &tHlts.L1_SingleMu16, &tHlts.b_L1_SingleMu16);
   t->SetBranchAddress("L1_SingleMu16_Prescl", &tHlts.L1_SingleMu16_Prescl, &tHlts.b_L1_SingleMu16_Prescl);
   t->SetBranchAddress("L1_SingleMu16er", &tHlts.L1_SingleMu16er, &tHlts.b_L1_SingleMu16er);
   t->SetBranchAddress("L1_SingleMu16er_Prescl", &tHlts.L1_SingleMu16er_Prescl, &tHlts.b_L1_SingleMu16er_Prescl);
   t->SetBranchAddress("L1_SingleMu20", &tHlts.L1_SingleMu20, &tHlts.b_L1_SingleMu20);
   t->SetBranchAddress("L1_SingleMu20_Prescl", &tHlts.L1_SingleMu20_Prescl, &tHlts.b_L1_SingleMu20_Prescl);
   t->SetBranchAddress("L1_SingleMu20er", &tHlts.L1_SingleMu20er, &tHlts.b_L1_SingleMu20er);
   t->SetBranchAddress("L1_SingleMu20er_Prescl", &tHlts.L1_SingleMu20er_Prescl, &tHlts.b_L1_SingleMu20er_Prescl);
   t->SetBranchAddress("L1_SingleMu25", &tHlts.L1_SingleMu25, &tHlts.b_L1_SingleMu25);
   t->SetBranchAddress("L1_SingleMu25_Prescl", &tHlts.L1_SingleMu25_Prescl, &tHlts.b_L1_SingleMu25_Prescl);
   t->SetBranchAddress("L1_SingleMu30", &tHlts.L1_SingleMu30, &tHlts.b_L1_SingleMu30);
   t->SetBranchAddress("L1_SingleMu30_Prescl", &tHlts.L1_SingleMu30_Prescl, &tHlts.b_L1_SingleMu30_Prescl);
   t->SetBranchAddress("L1_SingleMu5", &tHlts.L1_SingleMu5, &tHlts.b_L1_SingleMu5);
   t->SetBranchAddress("L1_SingleMu5_Prescl", &tHlts.L1_SingleMu5_Prescl, &tHlts.b_L1_SingleMu5_Prescl);
   t->SetBranchAddress("L1_SingleMu7", &tHlts.L1_SingleMu7, &tHlts.b_L1_SingleMu7);
   t->SetBranchAddress("L1_SingleMu7_Prescl", &tHlts.L1_SingleMu7_Prescl, &tHlts.b_L1_SingleMu7_Prescl);
   t->SetBranchAddress("L1_SingleMuBeamHalo", &tHlts.L1_SingleMuBeamHalo, &tHlts.b_L1_SingleMuBeamHalo);
   t->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &tHlts.L1_SingleMuBeamHalo_Prescl, &tHlts.b_L1_SingleMuBeamHalo_Prescl);
   t->SetBranchAddress("L1_SingleMuOpen", &tHlts.L1_SingleMuOpen, &tHlts.b_L1_SingleMuOpen);
   t->SetBranchAddress("L1_SingleMuOpen_Prescl", &tHlts.L1_SingleMuOpen_Prescl, &tHlts.b_L1_SingleMuOpen_Prescl);
   t->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &tHlts.L1_SingleMuOpen_NotBptxOR, &tHlts.b_L1_SingleMuOpen_NotBptxOR);
   t->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_Prescl", &tHlts.L1_SingleMuOpen_NotBptxOR_Prescl, &tHlts.b_L1_SingleMuOpen_NotBptxOR_Prescl);
   t->SetBranchAddress("L1_SingleJet128", &tHlts.L1_SingleJet128, &tHlts.b_L1_SingleJet128);
   t->SetBranchAddress("L1_SingleJet128_Prescl", &tHlts.L1_SingleJet128_Prescl, &tHlts.b_L1_SingleJet128_Prescl);
   t->SetBranchAddress("L1_SingleJet176", &tHlts.L1_SingleJet176, &tHlts.b_L1_SingleJet176);
   t->SetBranchAddress("L1_SingleJet176_Prescl", &tHlts.L1_SingleJet176_Prescl, &tHlts.b_L1_SingleJet176_Prescl);
   t->SetBranchAddress("L1_SingleJet200", &tHlts.L1_SingleJet200, &tHlts.b_L1_SingleJet200);
   t->SetBranchAddress("L1_SingleJet200_Prescl", &tHlts.L1_SingleJet200_Prescl, &tHlts.b_L1_SingleJet200_Prescl);
   t->SetBranchAddress("L1_SingleJet36", &tHlts.L1_SingleJet36, &tHlts.b_L1_SingleJet36);
   t->SetBranchAddress("L1_SingleJet36_Prescl", &tHlts.L1_SingleJet36_Prescl, &tHlts.b_L1_SingleJet36_Prescl);
   t->SetBranchAddress("L1_SingleJet52", &tHlts.L1_SingleJet52, &tHlts.b_L1_SingleJet52);
   t->SetBranchAddress("L1_SingleJet52_Prescl", &tHlts.L1_SingleJet52_Prescl, &tHlts.b_L1_SingleJet52_Prescl);
   t->SetBranchAddress("L1_SingleJet68", &tHlts.L1_SingleJet68, &tHlts.b_L1_SingleJet68);
   t->SetBranchAddress("L1_SingleJet68_Prescl", &tHlts.L1_SingleJet68_Prescl, &tHlts.b_L1_SingleJet68_Prescl);
   t->SetBranchAddress("L1_SingleJet92", &tHlts.L1_SingleJet92, &tHlts.b_L1_SingleJet92);
   t->SetBranchAddress("L1_SingleJet92_Prescl", &tHlts.L1_SingleJet92_Prescl, &tHlts.b_L1_SingleJet92_Prescl);
   t->SetBranchAddress("L1_SingleJetC20_NotBptxOR", &tHlts.L1_SingleJetC20_NotBptxOR, &tHlts.b_L1_SingleJetC20_NotBptxOR);
   t->SetBranchAddress("L1_SingleJetC20_NotBptxOR_Prescl", &tHlts.L1_SingleJetC20_NotBptxOR_Prescl, &tHlts.b_L1_SingleJetC20_NotBptxOR_Prescl);
   t->SetBranchAddress("L1_SingleJetC32_NotBptxOR", &tHlts.L1_SingleJetC32_NotBptxOR, &tHlts.b_L1_SingleJetC32_NotBptxOR);
   t->SetBranchAddress("L1_SingleJetC32_NotBptxOR_Prescl", &tHlts.L1_SingleJetC32_NotBptxOR_Prescl, &tHlts.b_L1_SingleJetC32_NotBptxOR_Prescl);
   t->SetBranchAddress("L1_TripleEG_14_10_8", &tHlts.L1_TripleEG_14_10_8, &tHlts.b_L1_TripleEG_14_10_8);
   t->SetBranchAddress("L1_TripleEG_14_10_8_Prescl", &tHlts.L1_TripleEG_14_10_8_Prescl, &tHlts.b_L1_TripleEG_14_10_8_Prescl);
   t->SetBranchAddress("L1_TripleMu0", &tHlts.L1_TripleMu0, &tHlts.b_L1_TripleMu0);
   t->SetBranchAddress("L1_TripleMu0_Prescl", &tHlts.L1_TripleMu0_Prescl, &tHlts.b_L1_TripleMu0_Prescl);
   t->SetBranchAddress("L1_TripleMu_5_5_3", &tHlts.L1_TripleMu_5_5_3, &tHlts.b_L1_TripleMu_5_5_3);
   t->SetBranchAddress("L1_TripleMu_5_5_3_Prescl", &tHlts.L1_TripleMu_5_5_3_Prescl, &tHlts.b_L1_TripleMu_5_5_3_Prescl);
   t->SetBranchAddress("L1_TripleJet_84_68_48_VBF", &tHlts.L1_TripleJet_84_68_48_VBF, &tHlts.b_L1_TripleJet_84_68_48_VBF);
   t->SetBranchAddress("L1_TripleJet_84_68_48_VBF_Prescl", &tHlts.L1_TripleJet_84_68_48_VBF_Prescl, &tHlts.b_L1_TripleJet_84_68_48_VBF_Prescl);
   t->SetBranchAddress("L1_TripleJet_92_76_64_VBF", &tHlts.L1_TripleJet_92_76_64_VBF, &tHlts.b_L1_TripleJet_92_76_64_VBF);
   t->SetBranchAddress("L1_TripleJet_92_76_64_VBF_Prescl", &tHlts.L1_TripleJet_92_76_64_VBF_Prescl, &tHlts.b_L1_TripleJet_92_76_64_VBF_Prescl);
   t->SetBranchAddress("L1_ZeroBias", &tHlts.L1_ZeroBias, &tHlts.b_L1_ZeroBias);
   t->SetBranchAddress("L1_ZeroBias_Prescl", &tHlts.L1_ZeroBias_Prescl, &tHlts.b_L1_ZeroBias_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0", &tHlts.L1Tech_BPTX_PreBPTX_v0, &tHlts.b_L1Tech_BPTX_PreBPTX_v0);
   t->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0_Prescl", &tHlts.L1Tech_BPTX_PreBPTX_v0_Prescl, &tHlts.b_L1Tech_BPTX_PreBPTX_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_minus.v0", &tHlts.L1Tech_BPTX_minus_v0, &tHlts.b_L1Tech_BPTX_minus_v0);
   t->SetBranchAddress("L1Tech_BPTX_minus.v0_Prescl", &tHlts.L1Tech_BPTX_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_minus_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0", &tHlts.L1Tech_BPTX_minus_AND_not_plus_v0, &tHlts.b_L1Tech_BPTX_minus_AND_not_plus_v0);
   t->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl", &tHlts.L1Tech_BPTX_minus_AND_not_plus_v0_Prescl, &tHlts.b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_plus.v0", &tHlts.L1Tech_BPTX_plus_v0, &tHlts.b_L1Tech_BPTX_plus_v0);
   t->SetBranchAddress("L1Tech_BPTX_plus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0", &tHlts.L1Tech_BPTX_plus_AND_NOT_minus_v0, &tHlts.b_L1Tech_BPTX_plus_AND_NOT_minus_v0);
   t->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &tHlts.L1Tech_BPTX_plus_AND_minus_v0, &tHlts.b_L1Tech_BPTX_plus_AND_minus_v0);
   t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_AND_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_AND_minus_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0", &tHlts.L1Tech_BPTX_plus_AND_minus_instance1_v0, &tHlts.b_L1Tech_BPTX_plus_AND_minus_instance1_v0);
   t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl", &tHlts.L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &tHlts.L1Tech_BPTX_plus_OR_minus_v0, &tHlts.b_L1Tech_BPTX_plus_OR_minus_v0);
   t->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_OR_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_OR_minus_v0_Prescl);
   t->SetBranchAddress("L1Tech_BPTX_quiet.v0", &tHlts.L1Tech_BPTX_quiet_v0, &tHlts.b_L1Tech_BPTX_quiet_v0);
   t->SetBranchAddress("L1Tech_BPTX_quiet.v0_Prescl", &tHlts.L1Tech_BPTX_quiet_v0_Prescl, &tHlts.b_L1Tech_BPTX_quiet_v0_Prescl);
   t->SetBranchAddress("L1Tech_CASTOR_Gap.v0", &tHlts.L1Tech_CASTOR_Gap_v0, &tHlts.b_L1Tech_CASTOR_Gap_v0);
   t->SetBranchAddress("L1Tech_CASTOR_Gap.v0_Prescl", &tHlts.L1Tech_CASTOR_Gap_v0_Prescl, &tHlts.b_L1Tech_CASTOR_Gap_v0_Prescl);
   t->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0", &tHlts.L1Tech_CASTOR_HaloMuon_v0, &tHlts.b_L1Tech_CASTOR_HaloMuon_v0);
   t->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_Prescl", &tHlts.L1Tech_CASTOR_HaloMuon_v0_Prescl, &tHlts.b_L1Tech_CASTOR_HaloMuon_v0_Prescl);
   t->SetBranchAddress("L1Tech_CASTOR_HighJet.v0", &tHlts.L1Tech_CASTOR_HighJet_v0, &tHlts.b_L1Tech_CASTOR_HighJet_v0);
   t->SetBranchAddress("L1Tech_CASTOR_HighJet.v0_Prescl", &tHlts.L1Tech_CASTOR_HighJet_v0_Prescl, &tHlts.b_L1Tech_CASTOR_HighJet_v0_Prescl);
   t->SetBranchAddress("L1Tech_CASTOR_MediumJet.v0", &tHlts.L1Tech_CASTOR_MediumJet_v0, &tHlts.b_L1Tech_CASTOR_MediumJet_v0);
   t->SetBranchAddress("L1Tech_CASTOR_MediumJet.v0_Prescl", &tHlts.L1Tech_CASTOR_MediumJet_v0_Prescl, &tHlts.b_L1Tech_CASTOR_MediumJet_v0_Prescl);
   t->SetBranchAddress("L1Tech_DT_GlobalOR.v0", &tHlts.L1Tech_DT_GlobalOR_v0, &tHlts.b_L1Tech_DT_GlobalOR_v0);
   t->SetBranchAddress("L1Tech_DT_GlobalOR.v0_Prescl", &tHlts.L1Tech_DT_GlobalOR_v0_Prescl, &tHlts.b_L1Tech_DT_GlobalOR_v0_Prescl);
   t->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0", &tHlts.L1Tech_HCAL_HBHE_totalOR_v0, &tHlts.b_L1Tech_HCAL_HBHE_totalOR_v0);
   t->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_Prescl", &tHlts.L1Tech_HCAL_HBHE_totalOR_v0_Prescl, &tHlts.b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl);
   t->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1", &tHlts.L1Tech_HCAL_HF_MMP_or_MPP_v1, &tHlts.b_L1Tech_HCAL_HF_MMP_or_MPP_v1);
   t->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1_Prescl", &tHlts.L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl, &tHlts.b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl);
   t->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2", &tHlts.L1Tech_HCAL_HF_coincidence_PM_v2, &tHlts.b_L1Tech_HCAL_HF_coincidence_PM_v2);
   t->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2_Prescl", &tHlts.L1Tech_HCAL_HF_coincidence_PM_v2_Prescl, &tHlts.b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl);
   t->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0", &tHlts.L1Tech_HCAL_HF_single_channel_v0, &tHlts.b_L1Tech_HCAL_HF_single_channel_v0);
   t->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0_Prescl", &tHlts.L1Tech_HCAL_HF_single_channel_v0_Prescl, &tHlts.b_L1Tech_HCAL_HF_single_channel_v0_Prescl);
   t->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0", &tHlts.L1Tech_HCAL_HO_totalOR_v0, &tHlts.b_L1Tech_HCAL_HO_totalOR_v0);
   t->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_Prescl", &tHlts.L1Tech_HCAL_HO_totalOR_v0_Prescl, &tHlts.b_L1Tech_HCAL_HO_totalOR_v0_Prescl);
   t->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RBplus1_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0);
   t->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl);
   t->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RBplus2_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0);
   t->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl);
   t->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_barrel_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_barrel_Cosmics_v0);
   t->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl);
   t->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_pointing_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_pointing_Cosmics_v0);
   t->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl);
   t->SetBranchAddress("L1Tech_TOTEM_0", &tHlts.L1Tech_TOTEM_0, &tHlts.b_L1Tech_TOTEM_0);
   t->SetBranchAddress("L1Tech_TOTEM_0_Prescl", &tHlts.L1Tech_TOTEM_0_Prescl, &tHlts.b_L1Tech_TOTEM_0_Prescl);
   t->SetBranchAddress("L1Tech_TOTEM_1", &tHlts.L1Tech_TOTEM_1, &tHlts.b_L1Tech_TOTEM_1);
   t->SetBranchAddress("L1Tech_TOTEM_1_Prescl", &tHlts.L1Tech_TOTEM_1_Prescl, &tHlts.b_L1Tech_TOTEM_1_Prescl);
   t->SetBranchAddress("L1Tech_TOTEM_2", &tHlts.L1Tech_TOTEM_2, &tHlts.b_L1Tech_TOTEM_2);
   t->SetBranchAddress("L1Tech_TOTEM_2_Prescl", &tHlts.L1Tech_TOTEM_2_Prescl, &tHlts.b_L1Tech_TOTEM_2_Prescl);
   t->SetBranchAddress("L1Tech_TOTEM_3", &tHlts.L1Tech_TOTEM_3, &tHlts.b_L1Tech_TOTEM_3);
   t->SetBranchAddress("L1Tech_TOTEM_3_Prescl", &tHlts.L1Tech_TOTEM_3_Prescl, &tHlts.b_L1Tech_TOTEM_3_Prescl);
   t->SetBranchAddress("L1Tech__TTU_RB0_Cosmics.v0", &tHlts.L1Tech__TTU_RB0_Cosmics_v0, &tHlts.b_L1Tech__TTU_RB0_Cosmics_v0);
   t->SetBranchAddress("L1Tech__TTU_RB0_Cosmics.v0_Prescl", &tHlts.L1Tech__TTU_RB0_Cosmics_v0_Prescl, &tHlts.b_L1Tech__TTU_RB0_Cosmics_v0_Prescl);
   if (doCheck) {
     if (t->GetMaximum("NL1IsolEm")>8) std::cout <<"FATAL ERROR: Arrary size of NL1IsolEm too small!!!  "<<t->GetMaximum("NL1IsolEm")<<std::endl;
     if (t->GetMaximum("NL1NIsolEm")>8) std::cout <<"FATAL ERROR: Arrary size of NL1NIsolEm too small!!!  "<<t->GetMaximum("NL1NIsolEm")<<std::endl;
     if (t->GetMaximum("NL1Mu")>8) std::cout <<"FATAL ERROR: Arrary size of NL1Mu too small!!!  "<<t->GetMaximum("NL1Mu")<<std::endl;
     if (t->GetMaximum("NL1CenJet")>8) std::cout <<"FATAL ERROR: Arrary size of NL1CenJet too small!!!  "<<t->GetMaximum("NL1CenJet")<<std::endl;
     if (t->GetMaximum("NL1ForJet")>8) std::cout <<"FATAL ERROR: Arrary size of NL1ForJet too small!!!  "<<t->GetMaximum("NL1ForJet")<<std::endl;
     if (t->GetMaximum("NL1Tau")>8) std::cout <<"FATAL ERROR: Arrary size of NL1Tau too small!!!  "<<t->GetMaximum("NL1Tau")<<std::endl;
   }
}

