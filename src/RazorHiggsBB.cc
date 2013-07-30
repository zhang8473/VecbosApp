// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "CommonTools/include/LumiReWeightingStandAlone.h"
#include "RazorHiggsBB.hh"

#define pileup_mc "2012S10MCPileUp.root"
#define pileup_data "pileup_190456-208686_8TeV_22Jan2013ReReco_Collisions12.root"
#define aJETPT 20.//another jet pT cut
#define JETPT 30.//jet pT cut
#define JETETA 2.5
#define HZ_MassMin 70.
#define CSVL 0.244
//#define CSVM 0.679
#define CSVM 0.5
#define CSVT 0.898
#define FIRSTCSVCUT CSVT
#define SECONDCSVCUT 0.5
#define MuonSelection muonPassTight
#define EleSelection electronPassWP95
#define MET2CUT 80*80
#define NUM_TREES 6
//#define MYDEBUG

RazorHiggsBB::RazorHiggsBB(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
}

RazorHiggsBB::RazorHiggsBB(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

void RazorHiggsBB::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorHiggsBB::SetWeight(double weight) {
  _weight = weight;
}

/*
bool SortBTag(const pair<TLorentzVector,Double_t> &j1, const pair<TLorentzVector,Double_t> &j2) {
  return ( j1.second>j2.second );
}

bool SortHZMass(const TLorentzVector &hz1, const TLorentzVector &hz2) {
    return ( hz1.M2()>hz2.M2() );
}
*/

void RazorHiggsBB::Loop(string outFileName, Long64_t start, Long64_t stop) {
  if(fChain == 0) return;
  reweight::LumiReWeighting LumiWeights( string(pileup_mc),string(pileup_data),"pileup","pileup" );

  /*  // prescaled RazorHiggsBB Triggers
  int HLT_R014_MR150;
  int HLT_R020_MR150;
  int HLT_R025_MR150;
  */
  //event weight
  Double_t Evt_Weight=1.;
  // hadronic razor triggers
  Bool_t HLT_R020_MR550,HLT_R025_MR450,HLT_R033_MR350,HLT_R038_MR250;
  
  //MET
  Double_t pfMET,DPhi_pfMET_trkMET,DPhi_pfMET_Jet;

  // general event info
  UChar_t Naj,Nal,N_totaljets,OrderInEvent;// number of other jets and leptons

  Double_t MAXBTAGAJet;// the max CVS of other jets
  //Double_t ColorPullAngleHZCands;
  Double_t BTAGHZ,jetareaHZ,jetpileupMVAHZ;//the merged HZ jet

  // prepare the output tree
  TTree* outTree[NUM_TREES];
  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  //TBranch* treeBranches[100];  UChar_t NBranches=0;
  char BranchTitle[100];
  outTree[0] = new TTree("CombinedJetsRazor", "CombinedJetsRazor");
  outTree[1] = new TTree("CombinedCSVLJetsRazor", "CombinedCSVLJetsRazor");
  outTree[2] = new TTree("CombinedCSVMJetsRazor", "CombinedCSVMJetsRazor");
  outTree[3] = new TTree("2BestCSVJetsRazor", "2BestCSVJetsRazor");
  outTree[4] = new TTree("2BJetsHZ", "2BJetsHZ");
  outTree[5] = new TTree("MergedJetsHZ", "MergedJetsHZ");

#define MakeBranch(Name,Var,Type) sprintf(BranchTitle,"%s/%s",Name,Type); \
  outTree[i]->Branch(Name, &Var,BranchTitle);

  for(int i=0; i<NUM_TREES; i++) {
    //event based selection

    MakeBranch("HLT_R020_MR550", HLT_R020_MR550, "O");
    MakeBranch("HLT_R025_MR450", HLT_R025_MR450, "O");
    MakeBranch("HLT_R033_MR350", HLT_R033_MR350, "O");
    MakeBranch("HLT_R038_MR250", HLT_R038_MR250, "O");
    
    if (_isData) {
      MakeBranch("run", runNumber, "I");
      MakeBranch("evNum", eventNumber, "l");
      MakeBranch("bx", bunchCrossing, "I");
      MakeBranch("ls", lumiBlock, "I");
      MakeBranch("orbit", orbitNumber, "I");
    }
    else {
      MakeBranch("evNum", eventNumber, "l");
      MakeBranch("nPU", nPU[1], "I");
    }

    MakeBranch("W", Evt_Weight, "D");
    MakeBranch("nPV", nPV, "I");

    //number of other jets Jets.size()-2
    MakeBranch("Naj", Naj, "b");
    //    MakeBranch("Nal", Nal, "b");

    //pfMET
    MakeBranch("pfMET", pfMET, "D");
    MakeBranch("pfMETphi", pfMETphi, "D");
    MakeBranch("DPhi_pfMET_trkMET", DPhi_pfMET_trkMET, "D");//azimuthal opening angle between the direction of ET measured with particle flow objects vs. charged particles (tkMET). This variable helps in reducing residual QCD background in the Z(νν)H channel arising from events where there is a mismatch in the missing energy measured by the calorime-ters vs. the tracker. It also has the advantage of being indepedent of ∆φ(pfMET, J).
    MakeBranch("DPhi_pfMET_Jet", DPhi_pfMET_Jet, "D");//azimuthal opening angle between the pfMET vector direction and the transverse momentum of the closest central jet in azimuth. Only jets satisfying pT > 30 GeV and |η | < 2.5 are considered. This variable helps in reducing residual QCD background in the Z(νν)H channel, where the source of the missing transverse energy is typically from fluctuations in the measured energy of a single jet.


    if (i<5) {
      //jet1
      MakeBranch("pTJet1", Jet1.pT, "D");
      MakeBranch("etaJet1", Jet1.eta, "D");
      MakeBranch("phiJet1", Jet1.phi, "D");
      MakeBranch("masssqJet1", Jet1.masssq, "D");
      
      //jet2
      MakeBranch("pTJet2", Jet2.pT, "D");
      MakeBranch("etaJet2", Jet2.eta, "D");
      MakeBranch("phiJet2", Jet2.phi, "D");
      MakeBranch("masssqJet2", Jet2.masssq, "D");
    }
    
    if (i==3||i==4) {//jet_others
      MakeBranch("BTAGJet1", Jet1.BTAG, "D");
      MakeBranch("pileupMVAJet1", Jet1.pileupMVA, "D");
      MakeBranch("areaJet1", Jet1.area, "D");
      MakeBranch("pullJet1", Jet1.pull, "D");
      MakeBranch("BTAGJet2", Jet2.BTAG, "D");
      MakeBranch("pileupMVAJet2", Jet2.pileupMVA, "D");
      MakeBranch("areaJet2", Jet2.area, "D");
      MakeBranch("pullJet2", Jet2.pull, "D");
    }

    if (i>0) MakeBranch("MAXBTAGAJet", MAXBTAGAJet, "D");

    if (i<5) {
      // MakeBranch("ColorPullAngle",ColorPullAngleHZCands,"D");
      MakeBranch("Deta_Jet1_Jet2",DEta_Jet1_Jet2,"D");
      MakeBranch("DR_Jet1_Jet2",DR_Jet1_Jet2,"D");
      //--R relating values and boost
      MakeBranch("MRsq", MRsq, "D");
      MakeBranch("RBeta", RBeta, "D");
      //--RStar relating values and boost
      MakeBranch("MRStarsq", MRStarsq, "D");
      MakeBranch("MTRsq", MTRsq, "D");
      MakeBranch("RStarsq", RStarsq, "D");
      MakeBranch("RStarBetaL", RStarBetaL, "D");
      MakeBranch("RStarBetaT", RStarBetaT, "D");      
    }

    //kinematic variables of the HZ candidate two daughters
    MakeBranch("pTHZ",   ptHZ, "D");
    MakeBranch("etaHZ",   etaHZ, "D");
    MakeBranch("phiHZ",   phiHZ, "D");
    MakeBranch("masssqHZ",   masssqHZ, "D");
    //topological positions of the HZ candidate two daughters
    MakeBranch("DPhi_H_Z",DPhi_pfMET_HZCands,"D");
    MakeBranch("DR_H_Z",DR_pfMET_HZCands,"D");
    if (i==5) {//merged HZ
	MakeBranch("BTAGHZ",   Jet1.BTAG, "D");
	MakeBranch("jetareaHZ",   Jet1.area, "D");
	MakeBranch("jetpileupMVAHZ",  Jet1.pileupMVA, "D");
	MakeBranch("jetpullHZ",  Jet1.pull, "D");
    }
    MakeBranch("OrderInEvent", OrderInEvent, "b");
  }

#define RENEWALL Jet1.pT=-9999.;Jet1.eta=-9999;Jet1.phi=-9999;Jet1.masssq=-9999;Jet1.BTAG=-9999;Jet1.pileupMVA=-9999;Jet1.area=-9999;Jet1.pull=-9999;Jet2.pT=-9999.;Jet2.eta=-9999;Jet2.phi=-9999;Jet2.masssq=-9999;Jet2.BTAG=-9999;Jet2.pileupMVA=-9999;Jet2.area=-9999;Jet2.pull=-9999;MAXBTAGAJet=-9999;DEta_Jet1_Jet2=-9999;DR_Jet1_Jet2=-9999;MRsq=-9999;RBeta=-9999;MRStarsq=-9999;MTRsq=-9999;RStarsq=-9999;RStarBetaL=-9999;RStarBetaT=-9999;ptHZ=-9999;etaHZ=-9999;phiHZ=-9999;masssqHZ=-9999;DPhi_pfMET_HZCands=-9999;DR_pfMET_HZCands=-9999;OrderInEvent=255

  //  Double_t _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // hadronic razor triggers
  std::vector<std::string> maskHLT_R020_MR550; maskHLT_R020_MR550.push_back("HLT_R020_MR550_v");
  std::vector<std::string> maskHLT_R025_MR450; maskHLT_R025_MR450.push_back("HLT_R025_MR450_v");
  std::vector<std::string> maskHLT_R033_MR350; maskHLT_R033_MR350.push_back("HLT_R033_MR350_v");
  std::vector<std::string> maskHLT_R038_MR250; maskHLT_R038_MR250.push_back("HLT_R038_MR250_v");

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Double_t sampleweight=1.;
  if (_isData) {
    Evt_Weight=1.;
    cout << "Number of entries = " << stop <<endl;
  }
  else {
    sampleweight=_weight/Double_t(stop);
    cout << "Number of entries = " << stop << "; xsec="<<_weight<< "; weight="<<sampleweight<<endl;
  }
  printf("Dummy");
  for (Long64_t jentry=start;  jentry<stop;jentry++) {//start event loop
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100 == 0||jentry>stop-20) printf("\r>>> Processing event # %ld/%ld",jentry+1,stop);

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      
      // hadronic razor triggers
      setRequiredTriggers(maskHLT_R020_MR550); reloadTriggerMask(true); HLT_R020_MR550 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R025_MR450); reloadTriggerMask(true); HLT_R025_MR450 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R033_MR350); reloadTriggerMask(true); HLT_R033_MR350 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R038_MR250); reloadTriggerMask(true); HLT_R038_MR250 = hasPassedHLT();
      
    }

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    // find highest-pT PV
    int iHighestPt = -99;
    Double_t HighestPt = -99999.;
    //if(nPV<1) continue;
    //for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    //if(ndofPV[iHighestPt] < 3) continue;
    //if(PVzPV[iHighestPt] > 25.) continue; 

    // HCAL FLAGS
    if(_isData && !eventPassHcalFilter()) continue;
    //ECALTPFilterFlag 
    /*
    if(_isData && METFlags << 0 == 0) continue;
    //drBoundary 
    if(_isData && METFlags << 1 == 0) continue;
    // drDead
    if(_isData && METFlags << 2 == 0) continue;
    // CSCHaloFilterFlag
    if(_isData && METFlags << 3 == 0) continue;
    // trackerFailureFilterFlag
    if(_isData && METFlags << 4 == 0) continue;
    // BE ECAL flag 
    if(_isData && METFlags << 5 == 0) continue;
    */
    // PFMET, PFMET cut
    _MET=TVector3(pxPFMet[0], pyPFMet[0], 0.);
    if ( _MET.Perp2()<MET2CUT ) continue;
    //lepton veto
    UShort_t NEle=0;
    for(UShort_t i=0; i<nEle; i++)
      if ( EleSelection(i) ) NEle++;
    UShort_t NMuon=0;
    for(UShort_t i=0; i<nMuon; i++)
      if ( MuonSelection(i) ) NMuon++;
    Nal=NEle+NMuon;
    if (Nal>0) continue;

    if (!_isData) Evt_Weight=sampleweight*LumiWeights.weight(nPU[1]);

    pfMET=_MET.Mag();
    pfMETphi=_MET.Phi();
    TVector3 trkMET(pxTCMet[0], pyTCMet[0], 0.);//Track-Corrected MET
    DPhi_pfMET_trkMET=OpeningAngle(pfMETphi,trkMET.Phi());
    DPhi_pfMET_Jet=99.;
    N_totaljets=0;
    // Jet selection 
    Bool_t goodPFevent=true;
    vector< pair<TLorentzVector,UShort_t> > Jets;
    for(int i=0; i< nAK5PFNoPUJet; i++) {//start jet loop
      Double_t pT2=pow(pxAK5PFNoPUJet[i],2.)+pow(pyAK5PFNoPUJet[i],2.);
      if( pT2<aJETPT*aJETPT ) continue;
      if( fabs(etaAK5PFNoPUJet[i])< 3.0 ) {//start eta 3.0 jets
	Double_t EU = uncorrEnergyAK5PFNoPUJet[i],
	  fnHAD = neutralHadronEnergyAK5PFNoPUJet[i]/EU,
	  fcHAD = chargedHadronEnergyAK5PFNoPUJet[i]/EU,
	  fnEM = neutralEmEnergyAK5PFNoPUJet[i]/EU,
	  fcEM = chargedEmEnergyAK5PFNoPUJet[i]/EU,
	  //fPhoton = photonEnergyAK5PFNoPUJet[i]/EU,
	  chargedMultiplicity = chargedHadronMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i],
	  neutralMultiplicity = neutralHadronMultiplicityAK5PFNoPUJet[i]+photonMultiplicityAK5PFNoPUJet[i];
	if ( fnHAD+fcHAD > 0.99 ) { goodPFevent = false; break; }//Only hadron scattering means a bad event
	if ( fnEM>0 || 
	     muonEnergyAK5PFNoPUJet[i] > 0.0 ||
	     ( fabs(etaAK5PFNoPUJet[i])  < JETETA && chargedEmEnergyAK5PFNoPUJet[i] > 0.0 ) ||
	     ( fabs(etaAK5PFNoPUJet[i])  < JETETA && chargedHadronEnergyAK5PFNoPUJet[i] > 0.0 )
	     ) {//Only forward electron/photon showering means a bad event
	  if (fnHAD < 0.99 && fnEM<0.99 && chargedMultiplicity+neutralMultiplicity>1 && fcHAD>0.0 && chargedMultiplicity>0 && fcEM<0.99 && fabs(etaAK5PFNoPUJet[i])< JETETA ) {
	    N_totaljets++;
	    if ( pT2>=JETPT*JETPT ) {
	      Float_t OpeningAnglewtMET=OpeningAngle(phiAK5PFNoPUJet[i],pfMETphi);
	      if (OpeningAnglewtMET<DPhi_pfMET_Jet) DPhi_pfMET_Jet=OpeningAnglewtMET;
	      Jets.push_back( make_pair(TLorentzVector(pxAK5PFNoPUJet[i], pyAK5PFNoPUJet[i], pzAK5PFNoPUJet[i], energyAK5PFNoPUJet[i]),i) );
	    }
	  }
	}	
	else {
	  goodPFevent = false;
	  clog<<"I found a bad PF event: Evt"<<eventNumber<<"--- fnEM:"<<fnEM<<";"<<endl;
	  break;
	}
      }//end eta 3.0 jets
    }//end jet loop
    if(goodPFevent == false || Jets.size()<2 ) continue;

    SortBTags(Jets);

    if (BTAG_Discriminator[Jets.front().second]<FIRSTCSVCUT) continue;
    if (BTAG_Discriminator[Jets[1].second]<SECONDCSVCUT) continue;

    //fill Branch 0: CombinedJetsRazor
    //RENEWALL;
    OrderInEvent=0;
    pair<TLorentzVector,TLorentzVector> CombinedJets = CombineJets(Jets,true);
    SetRazor(CombinedJets.first,CombinedJets.second);
    Naj=N_totaljets-Jets.size();
    outTree[0]->Fill();

    //fill Branch 1: CombinedCSVLJetsRazor
    //RENEWALL;
    vector< pair<TLorentzVector,UShort_t> >::iterator CSVL_Jet=Jets.begin();
    for (;CSVL_Jet!=Jets.end();CSVL_Jet++)
      if (BTAG_Discriminator[CSVL_Jet->second]<CSVL) break;
    if ( CSVL_Jet-Jets.begin()>1 ) {
      if ( CSVL_Jet!=Jets.end() ) MAXBTAGAJet=BTAG_Discriminator[CSVL_Jet->second];
      else MAXBTAGAJet=0;
#ifdef MYDEBUG
      cerr<<"njets="<<CSVL_Jet-Jets.begin()<<";CSVL_="<<BTAG_Discriminator[(CSVL_Jet-1)->second]<<endl;
#endif
      CombinedJets = CombineJets(vector< pair<TLorentzVector,UShort_t> >(Jets.begin(),CSVL_Jet),true);
      SetRazor(CombinedJets.first,CombinedJets.second);
      Naj=N_totaljets-(CSVL_Jet-Jets.begin());
      outTree[1]->Fill();
    }

    //fill Branch 2: CombinedCSVMJetsRazor
    //RENEWALL;
    vector< pair<TLorentzVector,UShort_t> >::iterator CSVM_Jet=Jets.begin();
    for (;CSVM_Jet!=Jets.end();CSVM_Jet++)
      if (BTAG_Discriminator[CSVM_Jet->second]<CSVM) break;
    if ( CSVM_Jet-Jets.begin()>1 ) {
      if ( CSVM_Jet!=Jets.end() ) MAXBTAGAJet=BTAG_Discriminator[CSVM_Jet->second];
      else MAXBTAGAJet=0;
#ifdef MYDEBUG
      cerr<<"njets="<<CSVM_Jet-Jets.begin()<<";CSVM_="<<BTAG_Discriminator[(CSVM_Jet-1)->second]<<endl;
#endif
      CombinedJets = CombineJets(vector< pair<TLorentzVector,UShort_t> >(Jets.begin(),CSVM_Jet),true);
      SetRazor(CombinedJets.first,CombinedJets.second);
      Naj=N_totaljets-(CSVM_Jet-Jets.begin());
      outTree[2]->Fill();
    }
    
    //fill Branch 3: 2BJetsRazor
    //RENEWALL;
    SetRazor(Jets[0].first,Jets[1].first);
    SetJets_Others(Jets[0].second,Jets[1].second);
    if (Jets.size()>2) MAXBTAGAJet=BTAG_Discriminator[Jets[2].second];
    else MAXBTAGAJet=0;
    Naj=N_totaljets-2;
    outTree[3]->Fill();
    
    //fill Branch 4: 2BJetsHZ
    //RENEWALL;
    OrderInEvent=0;
    for ( UChar_t i=0;i<Jets.size()-1;i++ )
      if ( BTAG_Discriminator[Jets[i].second]>FIRSTCSVCUT )
	for ( UChar_t j=i+1;j<Jets.size();j++ )
	  if ( BTAG_Discriminator[Jets[i].second]>SECONDCSVCUT )  {
	    TLorentzVector HZ = Jets[i].first+Jets[j].first;
	    if ( HZ.M()>HZ_MassMin ) {
	      SetRazor(Jets[i].first,Jets[j].first);
	      SetJets_Others(Jets[i].second,Jets[j].second);
	      ptHZ=HZ.Pt();
	      etaHZ=HZ.Eta();
	      phiHZ=HZ.Phi();
	      masssqHZ=HZ.M2();
	      if ( j+1<Jets.size() ) MAXBTAGAJet=BTAG_Discriminator[Jets[j+1].second];
	      else MAXBTAGAJet=0;
	      Naj=N_totaljets-2;
	      outTree[4]->Fill();
	      OrderInEvent++;
	    }
	  }
    
    //fill Branch 5: Merged HZ
    //RENEWALL;
    OrderInEvent=0;
    for ( UChar_t i=0;i<Jets.size();i++ ) 
      if (Jets[i].first.M() > HZ_MassMin) {
	ptHZ=Jets[i].first.Pt();
	etaHZ=Jets[i].first.Eta();
	phiHZ=Jets[i].first.Phi();
	masssqHZ=Jets[i].first.M2();
	SetJets_Others(Jets[i].second,0);
	DPhi_pfMET_HZCands=OpeningAngle(Jets[i].first.Phi(),_MET.Phi());
	DR_pfMET_HZCands=sqrt( pow(Jets[i].first.Eta()-_MET.Eta(),2)+DPhi_pfMET_HZCands*DPhi_pfMET_HZCands );
	if ( i+1<Jets.size() ) MAXBTAGAJet=BTAG_Discriminator[Jets[i+1].second];
	else MAXBTAGAJet=0;
	Naj=N_totaljets-1;
	outTree[5]->Fill();
	OrderInEvent++;
      }
  }//end the event loop
  printf("\nOutput to %s\n",outFileName.c_str());
  for(int i=0; i<NUM_TREES; i++) 
    outTree[i]->Write();
  file->Close();
}
  
pair<TLorentzVector,TLorentzVector> RazorHiggsBB::CombineJets(const vector< pair<TLorentzVector,UShort_t> > &myjets, const Bool_t SmallestSqMass){
  
  pair<TLorentzVector,TLorentzVector> CombinedJets;
  
  int N_comb = 1;
  for(int i = 0; i < myjets.size()-1; i++){
    N_comb *= 2;
  }
#ifdef MYDEBUG
  cerr<<endl<<myjets.size()<<" jets -- ";
#endif
  Double_t M_Extreme = SmallestSqMass?9999999999.0:0.;
  int j_count;
  for(int i = 1; i < N_comb; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i, count = 0;
    while ( itemp > 0) {
      if(itemp%2) j_temp1 += myjets[count].first;
      else j_temp2 += myjets[count].first;
#ifdef MYDEBUG
      cerr<<itemp%2;
#endif
      itemp/=2;
      count++;
    }
    for (vector< pair<TLorentzVector,UShort_t> >::const_iterator jet=myjets.begin()+count;jet!=myjets.end();jet++) {
      j_temp2+=jet->first;
#ifdef MYDEBUG
      cerr<<0;
#endif
    }
#ifdef MYDEBUG
    cerr<<":";
#endif
    Double_t M_temp = j_temp1.M2()+j_temp2.M2();
    if ( (SmallestSqMass && M_temp < M_Extreme) ||
	 (!SmallestSqMass && M_temp > M_Extreme)
	 ){
      M_Extreme = M_temp;
      CombinedJets.first = j_temp1;
      CombinedJets.second = j_temp2;
    }
  }

#ifdef MYDEBUG
  cerr<<endl;
#endif
  return CombinedJets;  
}
