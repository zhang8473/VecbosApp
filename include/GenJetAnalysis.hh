//-------------------------------------------------------
// Description:
//    Class for GenJet search analyses
// Authors:
//
//-------------------------------------------------------

/// The GenJetAnalysis class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef GenJetAnalysis_h
#define GenJetAnalysis_h

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <vector>

#include "VecbosBase.hh"
#include "Vecbos.hh"


class GenJetAnalysis : public Vecbos{
public:

  GenJetAnalysis(TTree *tree=0); /// Class Constructor
  virtual ~GenJetAnalysis();     /// Class Destructor
  /// The function to run on each events
  void Loop();
    
private:
  
  /// Create a set of histograms. Takes as input the name of the 
  /// directory were to write the histograms, also used to give 
  /// names to the histograms (to avoid memory problems if used
  /// for more than a set of histograms)
  vector<TH1D*> CreateHistos(string dirname); 

  /// Fill the histograms passes as input
  void FillHistos(vector<TH1D*> histos);

  double x0; ///< the X coordinate of the new vertex
  double y0; ///< the Y coordinate of the new vertex
  double z0; ///< the Z coordinate of the new vertex
};
#endif
