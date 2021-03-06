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
#define FIRSTJETPT 60.//first jet pT cut (only do jobs in 2BJetsHZ TTree)
#define JETETA 2.5
#define CSVL 0.244
#define CSVM 0.679
#define CSVT 0.898

/*Final Select Region
#define HZ_MassMin 60.
#define HZ_MassMax 120.
#define FIRSTCSVCUT CSVT
#define SECONDCSVCUT 0.5
#define EleSelection electronPassWP80
#define MET2CUT 100*100
#define MAXNLeptons 1
*/

//Including Control Region
#define HZ_MassMin 0.
#define HZ_MassMax 250.
#define FIRSTCSVCUT CSVL
#define SECONDCSVCUT 0
#define EleSelection electronPassWP80
#define MET2CUT 100*100
#define MAXNLeptons 1

#define NUM_TREES 6


const char *triggernames[]={"HLT_HT",
			    "HLT_RsqMR",
			    "HLT_PFNoPUHT",
			    "HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v",
			    "HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v",
			    "HLT_PFMET150_v",
			    "HLT_HT250_AlphaT0p55_v",
			    "HLT_PFNoPUHT350_PFMET100_v",
			    //"HLT_RsqMR40_Rsq0p04_v", prescaled trigger, do not use
			    "HLT_RsqMR55_Rsq0p09_MR150_v",
			    "HLT_RsqMR60_Rsq0p09_MR150_v",
			    "HLT_RsqMR65_Rsq0p09_MR150_v"};
const UChar_t Ntriggers=11;
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

Bool_t RazorHiggsBB::TauSelection(UShort_t itau) {
  return (isNonNullPFTau[itau]
	  &&abs(chargePFTau[itau])==1
	  &&thehpsTauDiscrByDecayModeFindingPFTau[itau]>0.5
	  &&pxPFTau[itau]*pxPFTau[itau]+pyPFTau[itau]*pyPFTau[itau]>1600
	  &&fabs(etaPFTau[itau])<2.1
	  &&theTauDiscrByLeadingTrackPtCutPFTau[itau]>0.5 //the pt cut is at 5Gev/c, discri has only 1 or 0
	  &&thehpsTauDiscrByTightMuonRejectionPFTau[itau]>0.5
	  &&thehpsTauDiscrByLooseElectronRejectionPFTau[itau]>0.5
	  &&thehpsTauDiscrByLooseCombinedIsolationDBSumPtCorrPFTau[itau]>0.5
	  )?true:false;//difference with AN-13-069: leading pt cut at 20GeV, one-prong only
}

Bool_t RazorHiggsBB::MuonSelection(UShort_t imu) {
  Utils AnalysisUtilities;
  int iTrack = trackIndexMuon[imu];
  Double_t PT = sqrt(pxMuon[imu] * pxMuon[imu] + pyMuon[imu] * pyMuon[imu]);
  return (AnalysisUtilities.muonIdVal(muonIdMuon[imu], bits::AllGlobalMuons)
	  &&AnalysisUtilities.muonIdVal(muonIdMuon[imu], bits::AllTrackerMuons)//AN-13-069(new): PFMuon
	  &&trackNormalizedChi2GlobalMuonTrack[iTrack]<10
	  &&numberOfValidPixelBarrelHitsTrack[iTrack]+numberOfValidPixelEndcapHitsTrack[iTrack]>0
	  &&numberOfValidStripTIBHitsTrack[iTrack]+numberOfValidStripTIDHitsTrack[iTrack]+numberOfValidStripTOBHitsTrack[iTrack]+numberOfValidStripTECHitsTrack[iTrack]>=6
	  &&numberOfValidMuonHitsGlobalMuonTrack[iTrack]>=1
	  &&numberOfMatchesMuon[imu]>=2
	  &&fabs(transvImpactParTrack[iTrack])<0.2
	  &&pfCandChargedIso04Muon[imu]+pfCandNeutralIso04Muon[imu]+pfCandPhotonIso04Muon[imu]<PT*0.12
	  &&fabs(etaMuon[imu])<2.4
	  &&PT>20
	  )?true:false;
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

  //event weight
  Double_t Evt_Weight=1.;
  // hadronic razor triggers
  Bool_t HasPassedHLT[Ntriggers+1];
  
  //MET
  Double_t pfMET,DPhi_pfMET_trkMET,DPhi_pfMET_Jet;

  // general event info
  UChar_t Naj,NEle,NMuon,NTau,N_totaljets,OrderInEvent;// number of other jets and leptons

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

  for (UChar_t i=0; i<NUM_TREES; i++) {
    //event based selection
    for (UChar_t j=0; j<Ntriggers;j++) {
      MakeBranch(triggernames[j], HasPassedHLT[j], "O");
    }
    MakeBranch("HLT_All_OR", HasPassedHLT[Ntriggers], "O");
    
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
      MakeBranch("W", Evt_Weight, "D");
    }

    MakeBranch("nPV", nPV, "I");

    //number of other jets Jets.size()-2
    MakeBranch("Naj", Naj, "b");
    MakeBranch("NEle", NEle, "b");
    MakeBranch("NMuon", NMuon, "b");
    MakeBranch("NTau", NTau, "b");

    //pfMET
    MakeBranch("pfMET", pfMET, "D");
    MakeBranch("pfMETsig", mEtSigMet[0], "D");
    MakeBranch("pfMETphi", pfMETphi, "D");
    MakeBranch("DPhi_pfMET_trkMET", DPhi_pfMET_trkMET, "D");//azimuthal opening angle between the direction of ET measured with particle flow objects vs. charged particles (tkMET). This variable helps in reducing residual QCD background in the Z(νν)H channel arising from events where there is a mismatch in the missing energy measured by the calorime-ters vs. the tracker. It also has the advantage of being indepedent of ∆φ(pfMET, J).
    MakeBranch("DPhi_pfMET_Jet", DPhi_pfMET_Jet, "D");//azimuthal opening angle between the pfMET vector direction and the transverse momentum of the closest central jet in azimuth. Only jets satisfying pT > 30 GeV and |η| < 2.5 are considered. This variable helps in reducing residual QCD background in the Z(νν)H channel, where the source of the missing transverse energy is typically from fluctuations in the measured energy of a single jet.


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
  std::vector<std::string> TriggerRequirements[Ntriggers+1];
  for (UChar_t i=0;i<Ntriggers;i++) {
    TriggerRequirements[i].push_back(std::string(triggernames[i]));
    TriggerRequirements[Ntriggers].push_back(std::string(triggernames[i]));
  }

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  if (_isData) cout << "Number of entries = " << stop <<endl;
  else cout << "Number of entries = " << stop << "; weight="<<_weight<<endl;

  printf("Dummy");
  Long64_t runNumber_prev_evt = 0;
  vector<int> requiredtriggerbits[Ntriggers+1];
  for (Long64_t jentry=start;  jentry<stop;jentry++) {//start event loop
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100 == 0||jentry>stop-20) printf("\r>>> Processing event # %ld/%ld",jentry+1,stop);

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    // hadronic razor triggers
    for (UChar_t i=0;i<Ntriggers+1;i++) {
      setRequiredTriggers(TriggerRequirements[i]);
      if ( (_isData&&runNumber!=runNumber_prev_evt)||jentry==start ) {
	//	cout<<"run "<<runNumber<<", reloading TriggerMask"<<endl;
	requiredtriggerbits[i]=reloadTriggerMask(true);
      }
      HasPassedHLT[i] = hasPassedHLT(requiredtriggerbits[i]);
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
    NEle=0;NMuon=0;NTau=0;
    for(UShort_t i=0; i<nEle; i++)
      if ( EleSelection(i) ) NEle++;
    for(UShort_t i=0; i<nMuon; i++)
      if ( MuonSelection(i) ) NMuon++;
    for(UShort_t i=0; i<nPFTau; i++)
      if ( TauSelection(i) ) NTau++;
    //printf("NEle=%d;NMu=%d;NTau=%d\n",NEle,NMuon,NTau);
    if (NEle+NMuon+NTau>MAXNLeptons) continue;

    if (!_isData) Evt_Weight=_weight*LumiWeights.weight(nPU[1]);

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
    for ( UChar_t i=0;i<Jets.size()-1;i++ )
      if ( BTAG_Discriminator[Jets[i].second]>FIRSTCSVCUT )
	for ( UChar_t j=i+1;j<Jets.size();j++ )
	  if ( BTAG_Discriminator[Jets[j].second]>SECONDCSVCUT &&
	       (Jets[i].first.Pt()>FIRSTJETPT||Jets[j].first.Pt()>FIRSTJETPT)
	       )  {
	    TLorentzVector HZ = Jets[i].first+Jets[j].first;
	    if ( HZ.M()>HZ_MassMin&&HZ.M()<HZ_MassMax ) {
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
      if (Jets[i].first.M() > HZ_MassMin && Jets[i].first.M()<HZ_MassMax) {
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
    runNumber_prev_evt=runNumber;
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
