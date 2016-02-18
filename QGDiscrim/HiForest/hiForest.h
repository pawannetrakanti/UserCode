#ifndef HIFOREST_H
#define HIFOREST_H

#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "SetupHltTree.h"
#include "SetupSkimTree.h"
#include "SetupEvtTree.h"
#include "SetupJetTree.h"
#include "SetupPFCandTree.h"
#include "SetupTrackTree.h"
#include "EventTree.h"
//#include "SetupTupleTree.h"

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TCut.h>

using namespace std;

typedef std::map<const char  *,int > JetMap;
typedef std::pair<const char *,int > JetPair;

// ==========================================================
// Main class which can be used to read the hiForest trees
//
// Author: Yen-Jie Lee
//
// ==========================================================

/* const int ndir=12; */
/* const char *calgo[ndir]= {"ak4Calo","ak3PF","ak4PF","ak5PF","akPu2Calo","akVs2Calo","akPu3Calo","akVs3Calo","akPu4Calo","akVs4Calo","akPu5Calo","akVs5Calo"}; */

const int ndir=6;
const char *calgo[ndir]= {"ak3PF","ak4PF","ak5PF","akPu3PF","akPu4PF","akPu5PF"};
//const char *calgo[ndir]= {"ak4Calo","ak3PF","ak4PF","ak5PF"};


class HiForest : public TObject
{

  using TObject::Draw;
  using TObject::GetName;

  public: 
  HiForest(const char *file, const char *name="forest", bool ispp = 0, bool ishireco=0, bool ismc = 0);
  virtual ~HiForest();

  // Utility functions
  Long64_t  GetEntry(int i);
  Long64_t  GetEntries();  						// Get the number of entries 
  void CheckTree(TTree *t,const char *title);			// Check the status of a tree
  void PrintStatus();						// Print the status of the hiForest

  Long64_t Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0){
    return tree->Draw(varexp,selection,option,nentries,firstentry);
  }
  const char *GetName(){
    return fname;
  }

  //! Jet Algorithms
  int GetNAlgo(); 
  const char *GetAlgoName(int i);
  Jets *GetJet(int i);
  Jets *GetJetByAlgo(const char *algo);
  void FillMap();
  void EmptyMap();
  void SelectJetAlgo (const char * algo[],const int len);
  void SelectBranches(const char *tname,const char * blst[],const int len);
  void SelectTupleBranches(std::vector<std::string> sblist);

  // TFile
  TFile *inf; 					// Input file 

  // Trees
  TTree *evtTree;                               // Event Tree
  TTree *hltTree;                               // Hlt Tree
  TTree *skimTree;                              // skim Tree
  TTree *pfTree;                                // PFCandidate Tree
  TTree *trackTree;                             // trackTree 
  TTree *JetTree[ndir];                         // Jet Trees
  TTree *tupleTree;                             // Tuple Tree

  TTree *tree;					// Pointer to the available tree, all trees in the forest are friended to each other

  // Branches

  // tree class for variables
  Hlts hlt; 
  Skims skim;
  Evts evt;
  PFCands pfcand;
  Tracks track;
  //Tuples *tuple;
  EventTree *tuple;

  //! Jet Trees
  Jets *mJets;
  JetMap JetContainer;
  JetMap::iterator it;

  // Booleans
  bool hasHltTree;
  bool hasSkimTree;
  bool hasPFCandTree;
  bool hasTrackTree;
  bool hasEvtTree;
  bool hasJetTree[ndir];
  bool hasTupleTree;

  bool pp;
  bool hireco;
  bool mc;

  Long64_t nEntries;
  Long64_t currentEvent;

 private:
  const char *fname;

  ClassDef(HiForest,9)
};
#endif
