//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorHiggsBB_h
#define RazorHiggsBB_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

#define BTAG_Discriminator combinedSecondaryVertexBJetTagsAK5PFNoPUJet

class RazorHiggsBB;
using namespace std;

class RazorHiggsBB : public Vecbos{
public:
  RazorHiggsBB() {}
  RazorHiggsBB(TTree *tree=0);// Class Constructor
  RazorHiggsBB(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorHiggsBB() {}     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, Long64_t start, Long64_t stop);
  void SetConditions(TTree* treeCond);
  Bool_t TauSelection(UShort_t itau);
  Bool_t MuonSelection(UShort_t iMuon);

  inline Double_t CalcMRsq(TLorentzVector ja, TLorentzVector jb){
    Double_t a0 = ja.P(), b0 = jb.P(),
      az = ja.Pz(), bz = jb.Pz();
    Double_t temp = a0*bz-b0*az;
    temp *= temp;
    temp /= (az-bz)*(az-bz)-(a0-b0)*(a0-b0);
    temp *= 4.;
    return temp;
  }
  inline Double_t CalcMRstarsq(TLorentzVector ja, TLorentzVector jb){
    Double_t a0 = ja.P(), b0 = jb.P(),
      az = ja.Pz(), bz = jb.Pz();
    TVector2 jaT(ja.Px(),ja.Py()), jbT(jb.Px(),jb.Py());
    return (a0+b0)*(a0+b0)-(az+bz)*(az+bz)-pow(jaT.Mod2()-jbT.Mod2(),2)/(jaT+jbT).Mod2();
  }
  inline Double_t CalcMTRsq(const TLorentzVector &ja, const TLorentzVector &jb, const TVector3 &met){
    return ( met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect()) )/2.;
  }
  inline Double_t CalcRStarBetaT(TLorentzVector ja, TLorentzVector jb){
    Double_t a0 = ja.P(), b0 = jb.P(),
      az = ja.Pz(), bz = jb.Pz();
    TVector2 jaT(ja.Px(),ja.Py()), jbT(jb.Px(),jb.Py());
    return fabs( jaT.Mod2()-jbT.Mod2() )/sqrt( (jaT+jbT).Mod2()*((a0+b0)*(a0+b0)-(az+bz)*(az+bz)) );
  }
  inline Double_t OpeningAngle(float phi1, float phi2) {//range from 0 to pi
    Double_t result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return fabs(result);
  }
private:

  pair<TLorentzVector,TLorentzVector> CombineJets(const vector< pair<TLorentzVector,UShort_t> > &myjets, const Bool_t IsSmallest);

  TVector3 _MET;
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  Double_t _weight;

  // HZ variables
  struct HZCand {
    Double_t Pt,Eta,Phi,OpenAnglewtMET,Mass;//the HZ candidate itself
    Double_t CSV1,CSV2,DEta,DR,ColorPullAngle;//the HZ daughters: if the two jets, all variables are meaningfull if one merged jet, CSV2,DEta,ColorPullAngle are not meaningfull;CSV1 is the jet CSV;DR is the area of the jet.
  };
  struct JetCand {
    Double_t pT,eta,phi,masssq,BTAG,pileupMVA,area,pull;
  }Jet1,Jet2;

  Double_t pfMETphi;
  Double_t MRsq,RBeta,MRStarsq,MTRsq,RStarsq,RStarBetaL,RStarGammaL,RStarBetaT;
  Double_t DEta_Jet1_Jet2,DR_Jet1_Jet2;//the topological info of the HZ candidate two daughters
  Double_t ptHZ,etaHZ,phiHZ,masssqHZ;//the HZ candidate
  Double_t DPhi_pfMET_HZCands,DR_pfMET_HZCands;//the topological info of the inv HZ and the visible ZH
  inline void SetRazor(TLorentzVector J1, TLorentzVector J2) {
    // the two jets
    if (J1.Pt() < J2.Pt() ) swap(J1,J2);
    Jet1.pT=J1.Pt();
    Jet1.eta=J1.Eta();
    Jet1.phi=J1.Phi();
    Jet1.masssq=J1.M2();
    Jet2.pT=J2.Pt();
    Jet2.eta=J2.Eta();
    Jet2.phi=J2.Phi();
    Jet2.masssq=J2.M2();
    DEta_Jet1_Jet2=Jet1.eta-Jet2.eta;
    DR_Jet1_Jet2=sqrt( pow(DEta_Jet1_Jet2,2)+pow(Jet1.phi-Jet2.phi,2) );
    //calculte the Razor variables
    //--R relating values and boost
    MRsq = CalcMRsq(J1, J2);
    RBeta = ( J1.E()-J2.E() )/( J1.Pz()-J2.Pz() );
    //--RStar relating values and boost
    MRStarsq = CalcMRstarsq(J1, J2);
    MTRsq = CalcMTRsq(J1, J2, _MET);
    RStarsq = MTRsq/MRStarsq;
     RStarBetaL =  ( J1.Pz()+J2.Pz() )/( J1.E()+J2.E() );
    if (RStarBetaL<=1) {
      Double_t RStarGammaL = sqrt( 1 - RStarBetaL*RStarBetaL );
    }
    else {cerr<<"BUG:TLorentzVec Pz>E"<<endl;}
    RStarBetaT = CalcRStarBetaT(J1, J2);
    DEta_Jet1_Jet2 = J1.Eta()-J2.Eta();
    DR_Jet1_Jet2=sqrt( pow(DEta_Jet1_Jet2,2)+pow(J1.Phi()-J2.Phi(),2) );
    TLorentzVector HZ = J1+J2;
    ptHZ=HZ.Pt();
    etaHZ=HZ.Eta();
    phiHZ=HZ.Phi();
    masssqHZ=HZ.M2();
    DPhi_pfMET_HZCands=OpeningAngle(HZ.Phi(),_MET.Phi());
    DR_pfMET_HZCands=sqrt( pow(HZ.Eta()-_MET.Eta(),2)+DPhi_pfMET_HZCands*DPhi_pfMET_HZCands );
  }

  inline void SetJets_Others(const UShort_t iJ1, const UShort_t iJ2) {
    Jet1.BTAG=combinedSecondaryVertexBJetTagsAK5PFNoPUJet[iJ1];Jet1.pileupMVA=jetIdMvaAK5PFNoPUJet[iJ1];
    Jet1.area=areaAK5PFNoPUJet[iJ1];Jet1.pull=pullAK5PFNoPUJet[iJ1];
    Jet2.BTAG=combinedSecondaryVertexBJetTagsAK5PFNoPUJet[iJ2];Jet2.pileupMVA=jetIdMvaAK5PFNoPUJet[iJ2];
    Jet2.area=areaAK5PFNoPUJet[iJ2];Jet2.pull=pullAK5PFNoPUJet[iJ2];
  }

  inline void SortBTags(vector< pair<TLorentzVector,UShort_t> > & jets) {
    UShort_t njets=jets.size();
    for (UShort_t i=0;i<njets-1;i++) 
      for (UShort_t j=i+1;j<njets;j++)
	if (BTAG_Discriminator[jets[i].second]<BTAG_Discriminator[jets[j].second])
	  swap(jets[i],jets[j]);
  }

};

#endif
