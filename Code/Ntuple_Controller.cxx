//Ntuple_Controller.cxx IMPLEMENTATION FILE
 

#include "Ntuple_Controller.h"
#include "Tools.h"
#include "PDG_Var.h"
#include "TF1.h"
#include "Parameters.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include <tuple>

// External code
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,cannotObtainHiggsMass(false)
  ,ObjEvent(-1)
  ,isInit(false)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("HTauTauTree");
  Logger(Logger::Verbose) << "Loading " << RootFiles.size() << " files" << std::endl;
  for(unsigned int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;  
  if(chain==0){
	Logger(Logger::Error) << "chain points to NULL" << std::endl;
  }
  Logger(Logger::Info) << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  Logger(Logger::Info) << "Ntuple Configured" << std::endl;


  // TFile *currentFile = chain->GetCurrentFile();
  // hLLRCounters = (TH1F*)currentFile->Get("HTauTauTree/Counters");
  // if(!hLLRCounters) std::cout<<"Counters histogram not found!"<<std::endl;
 


  // Fit setup 

  // Resolution uncertainty setup

  gRandom->SetSeed(1234);

  // Rochester muon momentum corrections

  rmcor = new rochcor2012(); // For systematics use rmcor = new rochcor2012(seed!=1234);

  // Set object correction flags to default values
  tauCorrection = "";
  muonCorrection = "";
  elecCorrection = "";
  jetCorrection = "";
}

///////////////////////////////////////////////////////////////////////
//
// Function: void InitEvent()
//
// Purpose: Initialize variables etc on event base
//
///////////////////////////////////////////////////////////////////////

void Ntuple_Controller::InitEvent(){
	Muon_corrected_p4.clear();
	//	Muon_corrected_p4.resize(NMuons());
	Muon_isCorrected = false;

	// after everything is initialized
	isInit = true;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
  isInit = false;
  InitEvent();
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  Logger(Logger::Verbose) << "Cleaning up" << std::endl;
  delete Ntp;
  delete rmcor;
  Logger(Logger::Verbose) << "Complete." << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
	Logger(Logger::Info) << "Starting D3PD cloning" << std::endl;
    newfile = new TFile(n+".root","recreate");
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  Logger(Logger::Info) << "Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  Logger(Logger::Warning) << "ThinTree not implemented." << std::endl;
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}


//Physics get Functions
Long64_t  Ntuple_Controller::GetMCID(){
	Long64_t  DataMCTypeFromTupel = Ntp->DataMC_Type_idx;
	//customize your event ID here 
	// if(DataMCTypeFromTupel==10230533 or DataMCTypeFromTupel==10130533 or DataMCTypeFromTupel==10330533 or DataMCTypeFromTupel==10430533) return 10230533;
	// if(DataMCTypeFromTupel==10410433 ) return DataMCTypeFromTupel;
	// move JAK Id information 3 digits to the left
	Long64_t  jakid = DataMCTypeFromTupel - (DataMCTypeFromTupel%100);
	jakid *= 1000;
	DataMCTypeFromTupel = jakid + (DataMCTypeFromTupel%100);
	if (DataMCTypeFromTupel == DataMCType::DY_ll_Signal && HistoC.hasID(DataMCType::DY_ll_Signal)) {
		for (unsigned int i = 0; i < NMCSignalParticles(); i++) {
			if (abs(MCSignalParticle_pdgid(i)) == PDGInfo::Z0) {
				if (fabs(MCSignalParticle_p4(i).M() - PDG_Var::Z_mass()) < 3 * PDG_Var::Z_width()) {
					return DataMCType::Signal;
				}
			}
		}
		return DataMCTypeFromTupel;
	}

	int dmcType = -999;

	// hack for Higgs mass splitting
	// Higgs mass is added to the MCId, such that the structure is JJJJJJAAABB (with JJJJJJ = JakID, AAA = mass, BB = DataMCType)

	if( (DataMCTypeFromTupel % 100) == DataMCType::H_tautau_ggF ||
	    (DataMCTypeFromTupel % 100) == DataMCType::H_tautau_VBF ||
	    (DataMCTypeFromTupel % 100) == DataMCType::H_tautau_WHZHTTH){
	  int mass = getSampleHiggsMass();
	  if (mass > 999)	Logger(Logger::Error) << "Read mass with more than 3 digits from sample: m = " << mass << std::endl;
	  if (mass > 0)		DataMCTypeFromTupel += mass*100;
	  // strip off JAK-Id from DataMCType
	  if (HistoC.hasID(DataMCTypeFromTupel % 100000)) {
	    dmcType = DataMCTypeFromTupel % 100000;
	  }
	}
	else {
	  // strip off JAK-Id from DataMCType
	  if (HistoC.hasID(DataMCTypeFromTupel % 100)) {
	    dmcType = DataMCTypeFromTupel % 100;
	  }
	}

	return dmcType;
}

// return DataMCType without mass information
int Ntuple_Controller::GetStrippedMCID(){
	return GetMCID() % 100;
}

// return path of input dataset (as given in the line InputNtuples in Input.txt)
TString Ntuple_Controller::GetInputNtuplePath(){
	Parameters Par; // assumes configured in Analysis.cxx
	TString dsPath;
	Par.GetString("InputNtuples:",dsPath);
	return dsPath;
}

// return name of input dataset
TString Ntuple_Controller::GetInputDatasetName(){
	TString dsPath = GetInputNtuplePath();
	return gSystem->BaseName( gSystem->DirName( gSystem->DirName(dsPath) ) );
}

//
TString Ntuple_Controller::GetInputPublishDataName(){
	TString dsPath = GetInputNtuplePath();
	return gSystem->BaseName( gSystem->DirName(dsPath) );
}

// determine Higgs mass (from Dataset name or fallback options)
int Ntuple_Controller::getSampleHiggsMass(){
	int mass = -999;

	// default method: analyze dataset name
	mass = readHiggsMassFromString( GetInputNtuplePath() );
	if (mass >= 0) return mass;

	// first fallback: analyze filename (only working when running on GRID)
	mass = readHiggsMassFromString( Get_File_Name() );
	if (mass >= 0) return mass;

	// second fallback: get Higgs mass from MC info
	Logger(Logger::Warning) << "Not able to obtain Higgs mass neither from dataset nor from file name."
	<< "\n\tMake sure the line InputNtuples is set correctly in your Input.txt. "
	<< "\n\tFor now, we will fall back to obtaining the Higgs mass from the generator information."
	<< "\n\tBe aware that SOME EVENTS WILL END UP IN THE WRONG HISTOGRAM!!!" << std::endl;
	return getHiggsSampleMassFromGenInfo();
}

int Ntuple_Controller::readHiggsMassFromString(TString input){
	// loop over possible masses
	for (int m = 100; m < 200; m = m+5){
		if ( input.Contains("M-" + TString::Itoa(m, 10) + "_") ) return m;
	}

	// mass not found
	return -999;
}

double Ntuple_Controller::getResonanceMassFromGenInfo(bool useZ0 /* = true*/, bool useHiggs0 /* = true*/, bool useW /* = true*/){
	for (unsigned int i = 0; i < NMCSignalParticles(); i++) {
		unsigned pdgid = abs(MCSignalParticle_pdgid(i));
		bool Z0 = useZ0 && (pdgid == PDGInfo::Z0);
		bool H0 = useHiggs0 && (pdgid == PDGInfo::Higgs0);
		bool W = useW && (pdgid == PDGInfo::W_plus);
		if (Z0 || H0 || W) {
			return MCSignalParticle_p4(i).M();
		}
	}
	return -999;
}

int Ntuple_Controller::getHiggsSampleMassFromGenInfo(){
	double resMass = getResonanceMassFromGenInfo(false, true, false);
	for (int m = 100; m < 200; m = m+5){
		if (fabs(resMass - m) < 2.5 ) {
			return m;
		}
	}
	return -999;
}

// TMatrixF Ntuple_Controller::Vtx_Cov(unsigned int i){
//   unsigned int dim=3;
//   TMatrixF M(dim,dim);
//   for(unsigned int j=0;j<dim;j++){
//     for(unsigned int k=0;k<=j;k++){
//       M[j][k]=Ntp->Vtx_Cov->at(i).at(j).at(k);
//       M[k][j]=Ntp->Vtx_Cov->at(i).at(j).at(k);
//     }
//   }
//   return M;
// }

// bool Ntuple_Controller::isVtxGood(unsigned int i){
//   if(0<=i && i<NVtx()){
//     if(Vtx_Track_idx(i).size()>4)return true;
//   }
//   return false;
// }

// bool Ntuple_Controller::isGoodVtx(unsigned int i){
// 	if(fabs(Vtx(i).z())>=24) return false;
// 	if(Vtx(i).Perp()>=2) return false;
// 	if(Vtx_ndof(i)<=4) return false;
// 	if(Vtx_isFake(i)) return false;
// 	return true;
// }

 // bool Ntuple_Controller::isGoodMuon(unsigned int i){
 //   //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
 //   //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
 //   //  isGoodMuon_nooverlapremoval(i) with
 //   //  ΔR(μ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
 //   if(isGoodMuon_nooverlapremoval(i)){
 //     unsigned int jet_idx=0;
 //     return !muonhasJetOverlap(i,jet_idx);
 //   }
 //   return false;
 // }

// bool Ntuple_Controller::muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx){
//   for(unsigned int j=0;j<NPFJets();j++){
//     if(isGoodJet_nooverlapremoval(j)){
//       if(Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))>0.2 && Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))<0.4){ jet_idx=j;return true;}
//     }
//   }
//   return false;
// }

// bool Ntuple_Controller::muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx){
//   for(unsigned int j=0;j<NPFJets();j++){
//     if(isGoodJet_nooverlapremoval(j)){
//       if(Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))<0.2){ jet_idx=j;return true;}
//     }
//   }
//   return false;
// }

 bool Ntuple_Controller::isTau(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   if(particleType(i)==2){
     if(Daughters_decayModeFindingOldDMs(i)>0.5){
	     return true;
     }
   }
   return false;
 }



 bool Ntuple_Controller::isLooseGoodTau(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   if(particleType(i)==2){

     //     std::cout<<

     if(Daughters_decayModeFindingOldDMs(i)>0.5)
       {
	 //       if(((tauID(i) & (1 << Bit_byLooseIsolationMVArun2v1DBoldDMwLT))==(1 << Bit_byLooseIsolationMVArun2v1DBoldDMwLT)))
	 {
	 if( ((tauID(i) & (1 << Bit_againstMuonLoose3))==(1 << Bit_againstMuonLoose3))){
	   if( ((tauID(i) & (1 << Bit_againstElectronVLooseMVA6))==(1 << Bit_againstElectronVLooseMVA6))){
	     return true;
	   }
	 }
       }
       }
   }
   return false;
 }



 // bool Ntuple_Controller::isLooseGoodTau(unsigned int i){
 //  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
 //   int tauIDmaskLoose(0);
 //   if(particleType(i)==2){
 //     tauIDmaskLoose |= (1<<Bit_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
 //     //     tauIDmaskLoose |= (1<<Bit_againstMuonLoose3);
 //     //     tauIDmaskLoose |= (1<<Bit_againstElectronVLooseMVA6); 
 //     std::cout<<" NC:  " << tauIDmaskLoose << "  dec   " << Daughters_decayModeFindingOldDMs(i)  <<"   "  <<tauID(i) <<std::endl;
 //     if(Daughters_decayModeFindingOldDMs(i) > 0.5){
 //       if(((tauID(i) & (1 <<tauIDmaskLoose ))==(1 << tauIDmaskLoose))){
     
 // 	 return true;
	 
 //       }
 //     }
 //   }
 //   return false;
 // }




 bool Ntuple_Controller::isMediumGoodTau(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   int tauIDmaskMedium(0);
   if(particleType(i)==2){
     tauIDmaskMedium|= (1<<Bit_againstMuonTight3);
     tauIDmaskMedium|= (1<<Bit_againstElectronMediumMVA6);
     if(Daughters_decayModeFindingOldDMs(i)>0.5){
       if((tauID(i) & tauIDmaskMedium) == tauIDmaskMedium){
     
	     return true;
	 }
       
     }
   }
   return false;
 }
  bool Ntuple_Controller::isTightGoodTau(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   int tauIDmaskMedium(0);
   if(particleType(i)==2){
     tauIDmaskMedium|= (1<<Bit_againstMuonTight3);
     tauIDmaskMedium|= (1<<Bit_againstElectronTightMVA6);
    

     if(Daughters_decayModeFindingOldDMs(i)>0.5){
       if((tauID(i) & tauIDmaskMedium) == tauIDmaskMedium){
	 return true;
       }
     }
   }
   return false;
 }
  


// bool Ntuple_Controller::isVTightIsolatedTau(int i){
//   int tauIDmaskMedium(0);
//   if(particleType(i)==2){
//     tauIDmaskMedium|= (1<<Bit_byVTightIsolationMVArun2v1DBoldDMwLT);
//     if((tauID(i) & tauIDmaskMedium) == tauIDmaskMedium){
//       return true;
//     }
//   }
//   return false;
// }



bool Ntuple_Controller::isIsolatedTau(int i, TString isotype){
  if(particleType(i)!=2){ std::cout<<"candidate is  not a tau  "<< i <<std::endl; return false;} 
  if(isotype.Contains("Loose")) return  CHECK_BIT(tauID(i),Bit_byLooseIsolationMVArun2v1DBoldDMwLT);
  if(isotype.Contains("Medium")) return  CHECK_BIT(tauID(i),Bit_byMediumIsolationMVArun2v1DBoldDMwLT);
  if(isotype.Contains("Tight")) return  CHECK_BIT(tauID(i),Bit_byTightIsolationMVArun2v1DBoldDMwLT);
  if(isotype.Contains("VTight")) return  CHECK_BIT(tauID(i),Bit_byVTightIsolationMVArun2v1DBoldDMwLT);
  return true;
}


bool Ntuple_Controller::tauBaselineSelection(int i, double cutPt, double cutEta, int aele, int amu){
  bool agEleVal(false), agMuVal(false), kin(false), dm(false), vertexS(false);
  // ag ele:
  if (aele== 0)      agEleVal = CHECK_BIT(tauID(i),Bit_againstElectronVLooseMVA6);
  else if ( aele== 1) agEleVal = CHECK_BIT(tauID(i),Bit_againstElectronLooseMVA6);
  else if ( aele== 2) agEleVal = CHECK_BIT(tauID(i),Bit_againstElectronMediumMVA6);
  else if ( aele== 3) agEleVal = CHECK_BIT(tauID(i),Bit_againstElectronTightMVA6);
  else if ( aele== 4) agEleVal = CHECK_BIT(tauID(i),Bit_againstElectronVTightMVA6);
  // ag mu:
  if ( amu== 0)      agMuVal = CHECK_BIT(tauID(i),Bit_againstMuonLoose3);
  else if ( amu== 1) agMuVal = CHECK_BIT(tauID(i),Bit_againstMuonTight3);

  kin     = (Daughters_P4(i).Pt() >cutPt && fabs(Daughters_P4(i).Eta())<cutEta );
  vertexS = (fabs(dz(i)) < 0.2); 
  dm      =(particleType(i)==2 && Daughters_decayModeFindingOldDMs(i) > 0.5);

  if(kin && agEleVal && agMuVal && dm && vertexS) return true;
  return false;
} 


bool Ntuple_Controller::muonBaselineSelection (int i, float ptMin, float etaMax, int muWP)
{
  bool kin(false),vertexS(false), idS(false), isMuon(false);
  if (muWP < 0 || muWP > 3)
    {
      cout << " ** NTuple_Controller::muBaseline: muWP must be between 0 and 3 --> using 0" << endl;
      muWP = 0;
    }
  isMuon = (particleType(i)==0);
  vertexS = (fabs(dxy(i)) < 0.045 && fabs(dz(i)) < 0.2);
  idS = CHECK_BIT(Daughters_muonID(i), muWP);
  kin     = (Daughters_P4(i).Pt() >ptMin && fabs(Daughters_P4(i).Eta())<etaMax );
  if(kin && vertexS && idS) return true;
  return false;
}


	   


std::vector<int>  Ntuple_Controller::SortTauHTauHPair (std::vector<int>  PairIndices)
{

  //  sorting according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#Baseline_tau_h_tau_h
  //  First prefer the pair with the most isolated candidate 1 (muon for μτh and eμ, electron for eτh and either τh for τhτh).
  //  If the isolation of candidate 1 is the same in both pairs, prefer the pair with the highest candidate 1 pt 
  //  (for cases of genuinely the same isolation value but different possible candidate 1).
  //  If the pt of candidate 1 in both pairs is the same (likely because it's the same object) then prefer the pair with the most isolated candidate 2 (tau for eτh, μτh, τhτh, electron for eμ).
  //  If the isolation of candidate 2 is the same, prefer the pair with highest candidate 2 pt (for cases of genuinely the same isolation value but different possible candidate 2). 
  //  std::vector<std::pair<double, double> >  vecpairs;

  // prepare a tuple with all the needed sorting info
  // pt1 - iso1 - idx1 - pt2 - iso2 - idx2 - idxoriginalPair


  std::vector<int> Sorted;
  vector<tauPair_t> vPairs;
  for(unsigned int ipair=0; ipair<PairIndices.size(); ipair++ )
    {
      TLorentzVector p4_1 =   Daughters_P4(indexDau1(ipair));
      TLorentzVector p4_2 =   Daughters_P4(indexDau2(ipair));
  
      float iso1 = Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(indexDau1(ipair));
      float iso2 = Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(indexDau2(ipair));
      // Daughters_byIsolationMVArun2v1DBoldDMwLTraw   - this is negative!
      // Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits  --  optional raw isolation
     // first one is highest pt
      tauPair_t pp;
      if (p4_1.Pt() > p4_2.Pt()) 
	{ 
	  pp.push_back(p4_1.Pt());  	pp.push_back(iso1); 	pp.push_back(indexDau1(ipair));
	  pp.push_back(p4_2.Pt()); 	pp.push_back(iso2); 	pp.push_back(indexDau2(ipair));	pp.push_back(PairIndices.at(ipair));
	  
	 //	pp = make_tuple(p4_1.Pt(), iso1, indexDau1(ipair), p4_2.Pt(), iso2,indexDau2(ipair) , ipair);
	}
       else { 	
        	pp.push_back(p4_2.Pt());	pp.push_back(iso2);	pp.push_back(indexDau2(ipair));
        	pp.push_back(p4_1.Pt());	pp.push_back(iso1);	pp.push_back(indexDau1(ipair));	pp.push_back(PairIndices.at(ipair));
      //  	//pp = make_tuple(p4_2.Pt(), iso2, indexDau2(ipair), p4_1.Pt(), iso1, indexDau1(ipair), ipair);
       }
      vPairs.push_back(pp);
    }
  std::vector<int> sortp;
  for(uint f = 0; f < vPairs.size(); ++f){
    int index(0);
    for(uint s = 0; s < vPairs.size(); ++s){
      if(s==f) continue;
      if(vPairs.at(f).at(1) >  vPairs.at(s).at(1)) {index = f;  continue;}
      else if(vPairs.at(f).at(0) > vPairs.at(s).at(0)){index = f;  continue;}
      else if(vPairs.at(f).at(4) > vPairs.at(s).at(4)){index = f;  continue;}
      else if(vPairs.at(f).at(3) > vPairs.at(s).at(3)){index = f;  continue;}
    }
    sortp.push_back(vPairs.at(index).at(6));
  }


   
  // now sort by iso, then pt criteria
  stable_sort(vPairs.begin(), vPairs.end(), pairSort);
  // PairIndices.clear();
   for(uint ipair = 0; ipair < vPairs.size(); ++ipair)
     {
       //       PairIndices.push_back( get<6> (vPairs.at(ipair)));
       Sorted.push_back(int(vPairs.at(ipair).at(6)));
     }
   return Sorted;
   //  return sortp;
}

bool Ntuple_Controller::pairSort (const tauPair_t& pA, const tauPair_t& pB)
{
  // first leg 1 iso
  float isoA = pA.at(1);////get<1> (pA);//
  float isoB = pB.at(1);////get<1> (pB);//
  if (isoA > isoB) return true; // NB: MVA iso ! Most iso --> highest MVA score
  else if (isoA < isoB) return false;

   // then leg 1 pt
   float ptA = pA.at(0);////get<0> (pA);//
   float ptB =pB.at(0);// //get<0> (pB);//
   if (ptA > ptB) return true;
   else if (ptA < ptB) return false;

   // then leg 2 iso
   isoA = pA.at(4);////get<4> (pA);//
   isoB = pB.at(4);////get<4> (pB);//
   if (isoA > isoB) return true;
   else if (isoA < isoB) return false;

   // then leg 2 pt
   ptA = pA.at(3);;////get<3> (pA);//
   ptB = pB.at(3);////get<3> (pB);//
   if (ptA > ptB) return true;
   else if (ptA < ptB) return false;

  // should be never here..
  return false;
}


bool Ntuple_Controller::pairSortRawIso (const tauPair_t& pA, const tauPair_t& pB)
{
  // first leg 1 iso
  float isoA = pA.at(1);//get<1> (pA);
  float isoB = pB.at(1);//get<1> (pB);//
  if (isoA < isoB) return true; // NB: raw iso ! Most iso --> lowest value
  else if (isoA > isoB) return false;

  // then leg 1 pt
  float ptA = pA.at(0);//get<0> (pA);//
  float ptB = pB.at(0);//get<0> (pB);//
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // then leg 2 iso
  isoA = pA.at(4);//get<4> (pA);//
  isoB = pB.at(4);//get<4> (pB);//
  if (isoA < isoB) return true;
  else if (isoA > isoB) return false;

  // then leg 2 pt
  ptA = pA.at(3);//get<3> (pA);//
  ptB = pB.at(3);//get<3> (pB);//
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  std::cout<<"  should be never here.. " <<std::endl;
  return false;
}








bool Ntuple_Controller::tauBaseline (int iDau, float ptMin, 
         float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits, 
         TString whatApply, bool debug)
{

  TLorentzVector p4 =Daughters_P4(iDau);

    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_dmfS  = false;
    bool byp_agEleS = false;
    bool byp_agMuS = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", 
    // "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All") && 
        !whatApply.Contains("SScharge") && 
        !whatApply.Contains("OScharge"))
    {
      byp_vertexS = byp_dmfS = byp_agEleS = byp_agMuS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_dmfS = false; 
      if (whatApply.Contains("againstEle"))  byp_agEleS = false; 
      if (whatApply.Contains("againstMu"))   byp_agMuS = false;       
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }


    if (againstEleWP < 0 || againstEleWP > 4) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstEleWP must be between 0 and 4 --> using 0" << endl;
        againstEleWP = 0;
    } 

    if (againstMuWP < 0 || againstMuWP > 1) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstMuWP must be between 0 and 1 --> using 0" << endl;
        againstMuWP = 0;
    }
    
    int agEleVal = 0;
    int agMuVal = 0;
    
    // ag ele:
    if (againstEleWP == 0)      agEleVal = CHECK_BIT(tauID(iDau),Bit_againstElectronVLooseMVA6);
    else if (againstEleWP == 1) agEleVal = CHECK_BIT(tauID(iDau),Bit_againstElectronLooseMVA6);
    else if (againstEleWP == 2) agEleVal = CHECK_BIT(tauID(iDau),Bit_againstElectronMediumMVA6);
    else if (againstEleWP == 3) agEleVal = CHECK_BIT(tauID(iDau),Bit_againstElectronTightMVA6);
    else if (againstEleWP == 4) agEleVal = CHECK_BIT(tauID(iDau),Bit_againstElectronVTightMVA6);

    // ag mu:
    if (againstMuWP == 0)      agMuVal = CHECK_BIT(tauID(iDau),Bit_againstMuonLoose3);
    else if (againstMuWP == 1) agMuVal = CHECK_BIT(tauID(iDau),Bit_againstMuonTight3);

    //bool dmfS = (tree->daughters_decayModeFindingOldDMs->at(iDau) == 1 || tree->daughters_decayModeFindingNewDMs->at(iDau) == 1) || byp_dmfS;
    bool dmfS = (Daughters_decayModeFindingOldDMs(iDau) == 1) || byp_dmfS;
    // bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    bool vertexS = (fabs(dz(iDau)) < 0.2) || byp_vertexS;
    bool agEleS = (agEleVal == 1) || byp_agEleS; 
    bool agMuS  = (agMuVal == 1) || byp_agMuS; 
    bool isoS = (Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(iDau) < isoRaw3Hits) || byp_isoS;
    if (whatApply.Contains ("InvertIzo")) isoS = !isoS ;

    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;

    bool totalS = (dmfS && vertexS && agEleS && agMuS && isoS && ptS && etaS);
    if (debug)
    {
      cout << "@ tau baseline" << endl;
      cout << " dmfS    "  << dmfS    << " skipped? " << byp_dmfS << endl;
      cout << " vertexS "  << vertexS << " skipped? " << byp_vertexS << endl;
      cout << " agEleS  "  << agEleS  << " skipped? " << byp_agEleS << endl;
      cout << " agMuS   "  << agMuS   << " skipped? " << byp_agMuS << endl;
      cout << " isoS    "  << isoS    << " skipped? " << byp_isoS << endl;
      cout << " ptS     "  << ptS     << " skipped? " << byp_ptS << endl;
      cout << " etaS    "  << etaS    << " skipped? " << byp_etaS << endl;
    }
    return totalS;    
}



bool Ntuple_Controller::muBaseline (int iDau, float ptMin, 
     float etaMax, float relIso, int muIDWP, TString whatApply, bool debug)
{

    TLorentzVector p4=Daughters_P4(iDau);
    int discr = Daughters_muonID(iDau);

    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_idS  = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", 
    // "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All") && 
        !whatApply.Contains("SScharge") && 
        !whatApply.Contains("OScharge"))
    {
      byp_vertexS = byp_idS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_idS = false; 
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }
      
    if (muIDWP < 0 || muIDWP > 3)
    {
        cout << " ** Ntuple_Controller::muBaseline: muIDWP must be between 0 and 3 --> using 0" << endl;
        muIDWP = 0;
    }

    bool vertexS = (fabs(dxy(iDau)) < 0.045 && fabs(dz(iDau)) < 0.2) || byp_vertexS;
    bool idS =  CHECK_BIT(discr, muIDWP) || byp_idS; // bit 0 is LOOSE id, bit 2 is MEDIUM mu id, bit 3 is TIGHT mu id
    bool isoS = (combreliso(iDau) < relIso) || byp_isoS;
    if (whatApply.Contains ("InvertIzo")) isoS = !isoS ;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    if (debug)
    {
      cout << "@ mu baseline" << endl;
      cout << " idS     "  << idS     << " skypped? " << byp_idS << endl;
      cout << " vertexS "  << vertexS << " skypped? " << byp_vertexS << endl;
      cout << " isoS    "  << isoS    << " skypped? " << byp_isoS << endl;
      cout << " ptS     "  << ptS     << " skypped? " << byp_ptS << endl;
      cout << " etaS    "  << etaS    << " skypped? " << byp_etaS << endl;
    }


    return totalS;
}


bool 
Ntuple_Controller::eleBaseline (int iDau, float ptMin, float etaMax, float relIso, int MVAIDflag, TString whatApply, bool debug)  //----------- to b e checked for a specific analysis
{ 

   
  TLorentzVector p4=Daughters_P4(iDau);
    // bypasser(s) and taker according to the string content
    bool byp_vertexS = false;
    bool byp_idS  = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", "againstMu", "Vertex", "SScharge"; separate various arguments with a semicolon
    if (!whatApply.Contains("All") && 
        !whatApply.Contains("SScharge") && 
        !whatApply.Contains("OScharge"))
    {
      byp_vertexS = byp_idS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_idS = false; 
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }

    bool vertexS = (fabs(dxy(iDau)) < 0.045 && fabs(dz(iDau)) < 0.2) || byp_vertexS;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;
    //bool idS = CHECK_BIT (Daughters_iseleCUT(iDau), 3) || byp_idS; // 3 is TIGHT ele id CUT BASED

    // bool idS = EleMVAID (tree->discriminator->at (iDau), tree->daughters_SCeta->at (iDau), p4.Pt (), MVAIDflag) || byp_idS ; // 2015/07/09 PG

    bool idS = Daughters_iseleBDT(iDau) || byp_idS; // use it in ntuples produced after 11 June 2015, contains tight WP bool  
    //bool idS = tightEleMVAID (tree->discriminator->at(iDau), TMath::Abs(p4.Eta())) || byp_idS; // APPROX! Using lepton eta and not super cluster eta, discriminator contains ele BDT  
    bool isoS = (combreliso(iDau) < relIso) || byp_isoS;
    if (whatApply.Contains ("InvertIzo")) isoS = !isoS ;
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);

    if (debug)
    {
      cout << "@ ele baseline" << endl;
      // cout << " debug: stored WP 80: " << tree->daughters_iseleWP80->at(iDau) << endl;
      // cout << " debug: stored RAW  : " << tree->discriminator->at (iDau) << " pt: " << p4.Pt() << " eta: " << p4.Eta() << endl;
      cout << " idS     "  << idS     << " skypped? " << byp_idS << endl;
      cout << " vertexS "  << vertexS << " skypped? " << byp_vertexS << endl;
      cout << " isoS    "  << isoS    << " skypped? " << byp_isoS << endl;
      cout << " ptS     "  << ptS     << " skypped? " << byp_ptS << endl;
      cout << " etaS    "  << etaS    << " skypped? " << byp_etaS << endl;
    }

    return totalS;
    
}

int Ntuple_Controller::getPairType (int type1, int type2)
{
    int nmu = 0;
    int nele = 0;
    int ntau = 0;
    
    if (isMuon (type1) )     nmu++;
    if (isElectron (type1) ) nele++;
    if (isTau (type1) )      ntau++;

    if (isMuon (type2) )     nmu++;
    if (isElectron (type2) ) nele++;
    if (isTau (type2) )      ntau++;

    if (nmu == 1 && nele == 0 && ntau == 1) return (int) MuHad;
    if (nmu == 0 && nele == 1 && ntau == 1) return (int) EHad;
    if (nmu == 0 && nele == 0 && ntau == 2) return (int) HadHad;
    if (nmu == 2 && nele == 0 && ntau == 0) return (int) MuMu;
    if (nmu == 0 && nele == 2 && ntau == 0) return (int) EE;
    if (nmu == 1 && nele == 1 && ntau == 0) return (int) EMu;
    
    return -1;
}

bool Ntuple_Controller::pairPassBaseline (int iPair, TString whatApply, bool debug)
{
    int dau1index = indexDau1(iPair);
    int dau2index = indexDau2(iPair);
    int type1 = particleType(dau1index);
    int type2 = particleType(dau2index);
    int pairType = getPairType (type1, type2);

    if (debug) cout << ".. checking baseline of pair " << iPair << " idx=(" << dau1index << "," << dau2index << ")" << endl;
        
    float dR = DeltaRDau(dau1index, dau2index);
    bool drMin = (dR > 0.1);
    if (!drMin && debug)
      cout << "failed dR min as dR=" << dR << endl;

    bool isOS = isOSCand(iPair);
    if (whatApply.Contains("OScharge") && !isOS) {
      if (debug) cout<<"check baseline: OSCharge failed"<<endl;
        return false; // do not even check the rest if requiring the charge
      }
    if (whatApply.Contains("SScharge") && isOS) {
      if (debug) cout<<"check baseline: SSCharge failed"<<endl;
        return false; // for the same sign selection at the moment full selection over SS pairs
      }

    // pairs are always ordered as: e mu | e tau | mu tau  (e < mu < tau)
    // if same type of particle, highest pt one is the first
    bool leg1=false;
    bool leg2=false;
    if (pairType == MuHad)
    {
        float tauIso = whatApply.Contains("TauRlxIzo") ? 7.0 : 3.0 ;
        leg1 = muBaseline (dau1index, 23., 2.1, 0.15, MuTight, whatApply, debug);
        leg2 = tauBaseline (dau2index, 20., 2.3, aeleVLoose, amuTight, tauIso, whatApply, debug);
    }

    if (pairType == EHad)
    {
        float tauIso = whatApply.Contains("TauRlxIzo") ? 7.0 : 3.0 ;
        leg1 = eleBaseline (dau1index, 27., 2.1, 0.1, EMVATight, whatApply, debug);
        leg2 = tauBaseline (dau2index, 20., 2.3, aeleTight, amuLoose, tauIso, whatApply, debug);
    }

    // ordered by pT and not by most isolated, but baseline asked in sync is the same...
    if (pairType == HadHad)
    {
        float tauIso = whatApply.Contains("TauRlxIzo") ? 7.0 : 2.0 ;
        leg1 = tauBaseline (dau1index, 40., 2.1, aeleVLoose, amuLoose, tauIso, whatApply, debug);
        leg2 = tauBaseline (dau2index, 40., 2.1, aeleVLoose, amuLoose, tauIso, whatApply, debug);
    }

    if (pairType == EMu)
    {
        // leg1 = eleBaseline (dau1index, 13., 0.15, EMVALoose, whatApply, debug);
        // leg2 = muBaseline ( dau2index, 9., 2.4, 0.15, MuTight, whatApply, debug);
    }
    
    // e e, mu mu are still preliminary (not from baseline)
    if (pairType == EE)
    {
      // leg1 = eleBaseline ( dau1index, 25., 0.15, EMVALoose, whatApply, debug);
      // leg2 = eleBaseline ( dau2index, 25., 0.15, EMVALoose, whatApply, debug);
    }
    
    if (pairType == MuMu)
    {
      leg1 = muBaseline ( dau1index, 10., 2.4, 0.15, MuTight, whatApply, debug);
      leg2 = muBaseline ( dau2index, 10., 2.4, 0.15, MuTight, whatApply, debug);
      bool leg1ER = muBaseline ( dau1index, 23., 2.1, 0.15, MuTight, whatApply, debug);
      bool leg2ER = muBaseline ( dau2index, 23., 2.1, 0.15, MuTight, whatApply, debug);
      
      //bool Is1in2p1 = leg1ER ;
      //bool Is2in2p1 = leg2ER ;
      return ( ((leg1ER && leg2) || (leg2ER && leg1)) && drMin );
    }
    
    bool result = (leg1 && leg2);
    if(!leg1 && debug){
      cout<<"check baseline: leg1 failed"<<endl;
    }
    if(!leg2 && debug){
      cout<<"check baseline: leg2 failed"<<endl;
    }
    
    if (leg1&&leg2&&debug)
      cout << "check baseline: leg1 AND leg2 ok" << endl;

    // bool drMin = (dR > 0.5);
    result = (result && drMin);

    return result;
}


float Ntuple_Controller::DeltaRDau(int dau1idx, int dau2idx)
{
  TLorentzVector v1, v2;
  v1 =Daughters_P4(dau1idx);
  v2 =Daughters_P4(dau2idx);
  return v1.DeltaR(v2);
}


 // bool Ntuple_Controller::muonBaselineSelection(int i){
 //  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
 //   if(particleType(i)==0){
 //     if(Daughters_P4(i).Pt()>25 && fabs(Daughters_P4(i).Eta())<2.1)
 //       if(fabs(dz(i))<0.2 &&  fabs(dxy(i))<0.045 )
 // 	 {
 // 	 return true;
 // 	 }
 //   }
 //   return false;
 // }

 bool Ntuple_Controller::electronBaselineSelection(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   if(particleType(i)==1){
     if(Daughters_P4(i).Pt()>17 && fabs(Daughters_P4(i).Eta())<2.4)
     	 return true;
   }
   return false;
 }

// electron
bool Ntuple_Controller::isElectron(int i){
  //  int bit(0);
  if(particleType(i)==1){
    {
      return true;
    }
  }
  return false;
}
bool Ntuple_Controller::ElectronVeto(int i){
  // ele.pt()         > 10                            and
  // fabs(ele.eta())  < 2.5                           and
  // fabs(dxy)        < 0.045                         and
  // fabs(dz)         < 0.2                           and
  // MVA ID 90% efficiency WP                         and
  // elec.passConversionVeto()            and 
  // elec.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS))) <=1            and  
  // iso              < 0.3
  if(particleType(i)==1){
    if(Daughters_P4(i).Pt()>10 && fabs(Daughters_P4(i).Eta())<2.5){
      if(fabs(dz(i))<0.2 &&  fabs(dxy(i))<0.045 ){
	if(Daughters_passConversionVeto(i)){
	  if(Daughters_eleMissingHits(i) <=1){
	    if(Daughters_iseleWP90(i)){
	      if(combreliso(i) < 0.3){
		return true;
	      }
	    }
	  }
	}
      }
    }
  }
  return false;
}

bool Ntuple_Controller::MuonVeto(int i){
  // muon.pt()        > 10                            and
  // fabs(muon.eta)   < 2.4                           and
  // fabs(dxy)        < 0.045                         and
  // fabs(dz)         < 0.2                           and
  // "Medium" (HIP safe) ID                                      and
  // iso(cone size 0.4)              < 0.3
  if(particleType(i)==0){
    if(Daughters_P4(i).Pt()>10 && fabs(Daughters_P4(i).Eta())<2.4){
      if(fabs(dz(i))<0.2 &&  fabs(dxy(i))<0.045 ){
	if(combreliso(i) < 0.3){
	  if(CHECK_BIT(Daughters_muonID(i), muIDWP::MuMedium) ){
	    return true;
	  }
	}
      }
    }
  }
  return false;
}


// muon
 bool Ntuple_Controller::isMuon(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   //   int bit(0);
   if(particleType(i)==0){
     //  if( ((Daughters_typeOfMuon(i) & (1<< 0)) == (1<< 0))  &&	 (((Daughters_typeOfMuon(i) & (1<< 0)) == (1<< 1))  ||  ((Daughters_typeOfMuon(i) & (1<< 0)) == (1<< 2)) )   )
       {
       return true;
     }
   }
   return false;
 }


//Loose muon
 bool Ntuple_Controller::isLooseGoodMuon(int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   //    int bit(0);
   if(particleType(i)==0){
     //  if( ((Daughters_typeOfMuon(i) & (1<< 0)) == (1<< 0))  &&	 (((Daughters_typeOfMuon(i) & (1<< 0)) == (1<< 1))  ||  ((Daughters_typeOfMuon(i) & (1<< 0)) == (1<< 2)) )   )
       {
       return true;
     }
   }
   return false;
 }




 bool Ntuple_Controller::isSoftGoodMuon( int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   if(particleType(i)==0){
     if(((Daughters_muonID(i) & (1<<Bit_MuonSoft)) == (1<<Bit_MuonSoft))){
       return true;
     }
   }
   return false;
 }
 bool Ntuple_Controller::isMediumGoodMuon( int i){
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   if(particleType(i)==0){
    if(((Daughters_muonID(i) & (1<<Bit_MuonMedium)) == (1<<Bit_MuonMedium))){
       return true;
     }
   }
   return false;
 }
 bool Ntuple_Controller::isTightGoodMuon( int i){
   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
   if(particleType(i)==0){
    if(((Daughters_muonID(i) & (1<<Bit_MuonTight)) == (1<<Bit_MuonTight))){
       return true;
     }
   }
   return false;
 }


// Global muon 	recoMu.isGlobalMuon()
// Normalized global-track χ2 	recoMu.globalTrack()->normalizedChi2() < 3
// Tracker-Standalone position match 	recoMu.combinedQuality().chi2LocalPosition < 12
// Kick finder 	recoMu.combinedQuality().trkKink < 20
// Segment compatibility 	muon::segmentCompatibility(recoMu) > 0.303




// void Ntuple_Controller::CorrectMuonP4(){
// 	if(isInit){
// 		for(unsigned int i=0;i<NMuons();i++){
// 			TLorentzVector mup4 = Muon_p4(i,"");
// 			int runopt = 0; // 0: no run-dependece
// 			float qter = 1.0; // 1.0: don't care about muon momentum uncertainty
// 			if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
// 				rmcor->momcor_mc(mup4,Muon_Charge(i),runopt,qter);
// 			}else{
// 				rmcor->momcor_data(mup4,Muon_Charge(i),runopt,qter);
// 			}
// 			Muon_corrected_p4.at(i) = mup4;
// 		}
// 		Muon_isCorrected = true;
// 	}else{
// 		Muon_isCorrected = false;
// 		Logger(Logger::Warning) << "No muon corrections applied" << std::endl;
// 	}
// }

/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Get muon four-vector
//
// Options:
//  - "roch": will correct the momentum in data and MC using Rochester corrections. For systematics look in constructor of this class.
//  - "scale": if you don't use momentum corrections, use this to estimate the systematics on the scale (only MC).
//             Add "down" in order to estimate it for downward variations.
//  - "res": if you don't use momentum corrections, use this to estimate systematics caused by momentum resolution (only MC)
//

// TLorentzVector Ntuple_Controller::Muon_p4(unsigned int i, TString corr){
// 	TLorentzVector vec = TLorentzVector(Ntp->Muon_p4->at(i).at(1),Ntp->Muon_p4->at(i).at(2),Ntp->Muon_p4->at(i).at(3),Ntp->Muon_p4->at(i).at(0));
// 	if (corr == "default") corr = muonCorrection;
// 	if(corr.Contains("roch")){
// 		if(!Muon_isCorrected){
// 			CorrectMuonP4();
// 		}
// 		if(Muon_isCorrected){
// 			vec = Muon_corrected_p4.at(i);
// 		}
// 		else{
// 			Logger(Logger::Warning) << "No muon corrections applied" << std::endl;
// 		}
// 	}
// 	if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
// 		if(corr.Contains("scale")){
// 			if(!corr.Contains("down")) vec.SetPerp(vec.Perp()*1.002);
// 			else vec.SetPerp(vec.Perp()*0.998);
// 		}else if(corr.Contains("res")){
// 			vec.SetPerp(gRandom->Gaus(vec.Perp(),vec.Perp()*0.006));
// 		}
// 		if(corr.Contains("met")){
// 			if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.002);
// 			else vec.SetPerp(vec.Perp() * 0.998);
// 		}
// 	}
// 	return vec;
// }

/////////////////////////////////////////////////////////////////////
//
// Official muon id code
//

// bool Ntuple_Controller::isTightMuon(unsigned int i){
// 	if(!Muon_isGlobalMuon(i)) return false;
// 	if(!Muon_isPFMuon(i)) return false;
// 	if(Muon_normChi2(i)>=10.) return false;
// 	if(Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
// 	if(Muon_numberOfMatchedStations(i)<=1) return false;
// 	if(Muon_numberofValidPixelHits(i)<=0) return false;
// 	if(Muon_trackerLayersWithMeasurement(i)<=5) return false;
// 	return true;
// }

// bool Ntuple_Controller::isTightMuon(unsigned int i, unsigned int j, TString corr){
// 	if(j<0 || j>=NVtx()) return false;
// 	if(!isTightMuon(i)) return false;
// 	if(dxy(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=0.2) return false;
// 	if(dz(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=0.5) return false;
// 	return true;
// }

// bool Ntuple_Controller::isLooseMuon(unsigned int i){
// 	if(!Muon_isPFMuon(i)) return false;
// 	if(!(Muon_isGlobalMuon(i) || Muon_isTrackerMuon(i))) return false;
// 	return true;
// }

// float Ntuple_Controller::Muon_RelIso(unsigned int i, TString corr){
// 	return (Muon_sumChargedHadronPt04(i)+std::max(0.,Muon_sumNeutralHadronEt04(i)+Muon_sumPhotonEt04(i)-0.5*Muon_sumPUPt04(i)))/Muon_p4(i,corr).Pt();
// }

/////////////////////////////////////////////////////////////////////

// bool Ntuple_Controller::isSelectedMuon(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr){
// 	if(j<0 || j>=NVtx()) return false;
// 	if(!isTightMuon(i)) return false;
// 	if(dxy(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=impact_xy) return false;
// 	if(dz(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=impact_z) return false;
// 	return true;
// }

/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Get electron four-vector
//
// Options:
//  - "res": use this to estimate the impact of the electron energy resolution on your result (only MC).
//           resolution in barrel (endcaps) for 2012: 1.6% (4.1%). uncertainty on resolution: 10%
//           RegEnergy corresponds to the electron energy after regression but before scale correction (data)
//           or smearing (MC)
//

// TLorentzVector Ntuple_Controller::Electron_p4(unsigned int i, TString corr){
// 	TLorentzVector vec = TLorentzVector(Ntp->Electron_p4->at(i).at(1),Ntp->Electron_p4->at(i).at(2),Ntp->Electron_p4->at(i).at(3),Ntp->Electron_p4->at(i).at(0));
// 	if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
// 		if (corr == "default") corr = elecCorrection;
// 		if(corr.Contains("scale") && Electron_RegEnergy(i)!=0){
// 			if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * (1+Electron_RegEnergyError(i)/Electron_RegEnergy(i)));
// 			else vec.SetPerp(vec.Perp() * (1-Electron_RegEnergyError(i)/Electron_RegEnergy(i)));
// 		}
// 		if(corr.Contains("res")){
// 			if(Electron_RegEnergy(i)>0){
// 				if(fabs(Electron_supercluster_eta(i))<1.479){
// 					if(corr.Contains("down")) vec.SetPerp(vec.Perp() * gRandom->Gaus(Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0144) / Electron_RegEnergy(i));
// 					else vec.SetPerp(gRandom->Gaus(vec.Perp() * Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0176) / Electron_RegEnergy(i));
// 				}
// 				else if(fabs(Electron_supercluster_eta(i))<2.5){
// 					if(corr.Contains("down")) vec.SetPerp(vec.Perp() * gRandom->Gaus(Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0369) / Electron_RegEnergy(i));
// 					else vec.SetPerp(gRandom->Gaus(vec.Perp() * Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0451) / Electron_RegEnergy(i));
// 				}
// 				else{
// 					Logger(Logger::Warning) << "Eta out of range: " << Electron_supercluster_eta(i) << ". Returning fourvector w/o smearing for resolution uncertainties." << std::endl;
// 				}
// 			}
// 			else{
// 				Logger(Logger::Warning) << "Energy <= 0: " << Electron_RegEnergy(i) << ". Returning fourvector w/o smearing for resolution uncertainties." << std::endl;
// 			}
// 		}
// 		if(corr.Contains("met")){
// 			if(fabs(Electron_supercluster_eta(i))<1.479){
// 				if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.006);
// 				else vec.SetPerp(vec.Perp() * 0.994);
// 			}
// 			else if(fabs(Electron_supercluster_eta(i))<2.5){
// 				if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.015);
// 				else vec.SetPerp(vec.Perp() * 0.985);
// 			}
// 			else{
// 				Logger(Logger::Warning) << "Eta out of range: " << Electron_supercluster_eta(i) << ". Returning fourvector w/o scale corrections for MET uncertainties." << std::endl;
// 			}
// 		}
// 	}
// 	return vec;
// }

/////////////////////////////////////////////////////////////////////
//
// Official electron id code
//

// bool Ntuple_Controller::isTrigPreselElectron(unsigned int i){
// 	if(fabs(Electron_supercluster_eta(i))>2.5) return false;
// 	if(Electron_numberOfMissedHits(i)>0) return false;
// 	if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
// 	if(Electron_Gsf_dr03EcalRecHitSumE(i)/Electron_p4(i).Pt()>0.2) return false;
// 	if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
// 	if(fabs(Electron_supercluster_eta(i))<1.479){
// 		if(Electron_sigmaIetaIeta(i)>0.014) return false;
// 		if(Electron_hadronicOverEm(i)>0.15) return false;
// 	}else{
// 		if(Electron_sigmaIetaIeta(i)>0.035) return false;
// 		if(Electron_hadronicOverEm(i)>0.1) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isMVATrigElectron(unsigned int i, TString corr){
// 	// !!! make sure to also apply Electron_RelIso<0.15 in your analysis !!!
// 	double mvapt = Electron_p4(i,corr).Pt();
// 	double mvaeta = fabs(Electron_supercluster_eta(i));
// 	if(mvapt<10.) return false;
// 	if(mvaeta>2.5) return false;
// 	if(Electron_numberOfMissedHits(i)>0) return false;
// 	if(Electron_HasMatchedConversions(i)) return false;
// 	if(!isTrigPreselElectron(i)) return false;
// 	if(mvapt>10. && mvapt<20.){
// 		if(mvaeta<0.8 && Electron_MVA_Trig_discriminator(i)<=0.00) return false;
// 		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_Trig_discriminator(i)<=0.10) return false;
// 		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_Trig_discriminator(i)<=0.62) return false;
// 	}
// 	if(mvapt>=20.){
// 		if(mvaeta<0.8 && Electron_MVA_Trig_discriminator(i)<=0.94) return false;
// 		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_Trig_discriminator(i)<=0.85) return false;
// 		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_Trig_discriminator(i)<=0.92) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isTrigNoIPPreselElectron(unsigned int i){
// 	if(fabs(Electron_supercluster_eta(i))>2.5) return false;
// 	if(Electron_numberOfMissedHits(i)>0) return false;
// 	if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
// 	if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
// 	if(fabs(Electron_supercluster_eta(i))<1.479){
// 		if(Electron_sigmaIetaIeta(i)>0.01) return false;
// 		if(Electron_hadronicOverEm(i)>0.12) return false;
// 		if(fabs(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.007) return false;
// 		if(fabs(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.15) return false;
// 	}else{
// 		if(Electron_sigmaIetaIeta(i)>0.03) return false;
// 		if(Electron_hadronicOverEm(i)>0.1) return false;
// 		if(fabs(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.009) return false;
// 		if(fabs(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.1) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isMVATrigNoIPElectron(unsigned int i, TString corr){
// 	// at present there are no recommendations on the isolation
// 	double mvapt = Electron_p4(i, corr).Pt();
// 	double mvaeta = fabs(Electron_supercluster_eta(i));
// 	if(mvaeta>2.5) return false;
// 	if(Electron_HasMatchedConversions(i)) return false;
// 	if(Electron_numberOfMissedHits(i)>0) return false;
// 	if(!isTrigNoIPPreselElectron(i)) return false;
// 	if(mvapt<20){
// 		if(mvaeta<0.8 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.5375) return false;
// 		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.375) return false;
// 		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.025) return false;
// 	}
// 	if(mvapt>=20){
// 		if(mvaeta<0.8 && Electron_MVA_TrigNoIP_discriminator(i)<=0.325) return false;
// 		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
// 		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isMVANonTrigElectron(unsigned int i, unsigned int j, TString corr){
// 	// !!! make sure to also apply Electron_RelIso<0.4 in your analysis !!!
// 	double mvapt = Electron_p4(i,corr).Pt();
// 	double mvaeta = fabs(Electron_supercluster_eta(i));
// 	if(mvapt<7.) return false;
// 	if(mvaeta>2.5) return false;
// 	if(Electron_numberOfMissedHits(i)>1) return false;
// 	if(vertexSignificance(Electron_Poca(i),j)>=4) return false;
// 	if(mvapt>7. && mvapt<10.){
// 		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.47) return false;
// 		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.004) return false;
// 		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.295) return false;
// 	}
// 	if(mvapt>=10.){
// 		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=-0.34) return false;
// 		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=-0.65) return false;
// 		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.6) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isTightElectron(unsigned int i, TString corr){
// 	if(Electron_HasMatchedConversions(i)) return false;
// 	if(Electron_numberOfMissedHits(i)>0) return false;
// 	if(Electron_RelIso04(i,corr)>=0.1) return false;
// 	if(fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))>=0.05) return false;
// 	if(fabs(Electron_supercluster_eta(i))<=1.479){
// 		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.004) return false;
// 		if(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.03) return false;
// 		if(Electron_sigmaIetaIeta(i)>=0.01) return false;
// 		if(Electron_hadronicOverEm(i)>=0.12) return false;
// 	}
// 	if(fabs(Electron_supercluster_eta(i))>1.479 && fabs(Electron_supercluster_eta(i))<2.5){
// 		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.005) return false;
// 		if(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.02) return false;
// 		if(Electron_sigmaIetaIeta(i)>=0.03) return false;
// 		if(Electron_hadronicOverEm(i)>=0.10) return false;
// 		if(Electron_p4(i).Pt()<20 && Electron_RelIso04(i)>=0.07) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isTightElectron(unsigned int i, unsigned int j, TString corr){
// 	if(j<0 || j>=NVtx()) return false;
// 	if(!isTightElectron(i,corr)) return false;
// 	if(dxy(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=0.02) return false;
// 	if(dz(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=0.1) return false;
// 	return true;
// }

// float Ntuple_Controller::Electron_RelIso03(unsigned int i, TString corr){
// 	return (Electron_chargedHadronIso(i)+std::max((double)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff_R03(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
// }

// double Ntuple_Controller::Electron_RelIsoDep03(unsigned int i, TString corr){
// 	return (Electron_isoDeposits_chargedHadronIso03(i)+std::max((double)0.,Electron_isoDeposits_neutralHadronIso03(i)+Electron_isoDeposits_photonIso03(i)-RhoIsolationAllInputTags()*Electron_Aeff_R03(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
// }

// float Ntuple_Controller::Electron_RelIso04(unsigned int i, TString corr){
// 	return (Electron_chargedHadronIso(i)+std::max((double)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff_R04(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
// }

// double Ntuple_Controller::Electron_RelIsoDep04(unsigned int i, TString corr){
// 	return (Electron_isoDeposits_chargedHadronIso04(i)+std::max((double)0.,Electron_isoDeposits_neutralHadronIso04(i)+Electron_isoDeposits_photonIso04(i)-RhoIsolationAllInputTags()*Electron_Aeff_R04(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
// }

// double Ntuple_Controller::Electron_Aeff_R04(double Eta){
// 	double eta=fabs(Eta);
// 	if(eta>=0. && eta<1.) return 0.208;
// 	else if(eta>=1. && eta<1.479) return 0.209;
// 	else if(eta>=1.479 && eta<2.) return 0.115;
// 	else if(eta>=2. && eta<2.2) return 0.143;
// 	else if(eta>=2.2 && eta<2.3) return 0.183;
// 	else if(eta>=2.3 && eta<2.4) return 0.194;
// 	else if(eta>=2.4) return 0.261;
// 	else {Logger(Logger::Error) << "Electron eta out of range: " << Eta << std::endl; return -1;}
// }

// double Ntuple_Controller::Electron_Aeff_R03(double Eta){
// 	double eta=fabs(Eta);
// 	if(eta>=0. && eta<1.) return 0.13;
// 	else if(eta>=1. && eta<1.479) return 0.14;
// 	else if(eta>=1.479 && eta<2.) return 0.07;
// 	else if(eta>=2. && eta<2.2) return 0.09;
// 	else if(eta>=2.2 && eta<2.3) return 0.11;
// 	else if(eta>=2.3 && eta<2.4) return 0.11;
// 	else if(eta>=2.4) return 0.14;
// 	else {Logger(Logger::Error) << "Electron eta out of range: " << Eta << std::endl; return -1;}
// }

// /////////////////////////////////////////////////////////////////////

// bool Ntuple_Controller::isSelectedElectron(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr){
// 	double mvapt = Electron_p4(i,corr).Pt();
// 	double mvaeta = fabs(Electron_supercluster_eta(i));
// 	if(Electron_numberOfMissedHits(i)>0) return false;
// 	if(Electron_HasMatchedConversions(i)) return false;
// 	if(dxy(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=impact_xy) return false;
// 	if(dz(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=impact_z) return false;
// 	if(mvapt<20.){
// 		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.925) return false;
// 		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.915) return false;
// 		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.965) return false;
// 	}
// 	if(mvapt>=20.){
// 		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.905) return false;
// 		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.955) return false;
// 		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.975) return false;
// 	}
// 	return true;
// }

// bool Ntuple_Controller::isGoodJet(unsigned int i){
//   //  Top Dilepton Jet selection with pt 15GeV
//   //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel 
//   //  isGoodJet_nooverlapremoval(i) with:
//   //  deltaR jet-electron cleaning < 0.4 (2 selected lepton only) < 0.3 (e+jets only) 0.3
//   //  deltaR jet-muon cleaning < 0.4(2 selected lepton only) < 0.3 (mu+jets only), 0.1 for PF and JET 0.3
//   if(isGoodJet_nooverlapremoval(i)){
//     unsigned int muon_idx;
//     return !jethasMuonOverlap(i,muon_idx);
//   }
//   return false;
// }

// bool Ntuple_Controller::jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx){
//   for(unsigned int j=0;j<NMuons();j++){
//     if(isGoodMuon_nooverlapremoval(j) && Muon_RelIso(j)<0.2){
//       if(Tools::dr(Muon_p4(j),PFJet_p4(jet_idx))<0.4){ muon_idx=j;return true;}
//     }
//   }
//   return false;
// }


// bool Ntuple_Controller::isGoodJet_nooverlapremoval(unsigned int i){
//   //  Top Dilepton Jet selection with pt 15GeV
//   //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel
//   //  Jet ID defined in isJetID(i)               
//   //  cut         dilepton    l+jet
//   //  corrected pT 30 GeV 30 GeV    
//   //  residual correction (data) applied applied
//   //  abs(eta) < 2.5 < 2.4                      
//   //  jet ID applied applied     
//   if(isJetID(i)){
//     if(PFJet_p4(i).Pt()>15.0){
//       if(fabs(PFJet_p4(i).Eta())<2.4){
// 	return true;
//       }
//     }
//   }
//   return false;
// }

// bool Ntuple_Controller::isJetID(unsigned int i, TString corr){
//   //  Top Dilepton Jet selection with pt and iso matching the muon and tau.
//   //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
//   //  Jet ID :
//   //  corrected jet pt > 10 GeV, jet eta < 5.2
//   //  number of constituents>1 (patJet->numberOfDaughters())
//   //  NHF<0.99 (( patJet->neutralHadronEnergy() + patJet->HFHadronEnergy() ) / patJet->energy())
//   //  NEF<0.99 (patJet->neutralEmEnergyFraction())
//   //  if |η|<2.4, CEF<0.99 (patJet->chargedEmEnergyFraction())
//   //  if |η|<2.4, CHF>0 (patJet->chargedHadronEnergyFraction())
//   //  if |η|<2.4, NCH>0 (patJet->chargedMultiplicity()) 
//   /////////////////////////////////////////////////////////////////////////
//   // apply jet ID
//   if(PFJet_p4(i,corr).Pt()<=10.) return false;
//   if(fabs(PFJet_p4(i,corr).Eta())>=5.2) return false;
//   if(PFJet_numberOfDaughters(i)<=1) return false;
//   if((PFJet_neutralHadronEnergy(i)+PFJet_HFHadronEnergy(i))/PFJet_p4(i,corr).E()>=0.99) return false;
//   if(PFJet_neutralEmEnergyFraction(i)>=0.99) return false;
//   if(fabs(PFJet_p4(i,corr).Eta())<2.4){
// 	  if(PFJet_chargedEmEnergyFraction(i)>=0.99) return false;
// 	  if(PFJet_chargedHadronEnergyFraction(i)<=0.) return false;
// 	  if(PFJet_chargedMultiplicity(i)<=0) return false;
//   }
//   return true;
// }

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECL2ResidualTimeStability#2012Rereco
// double Ntuple_Controller::rundependentJetPtCorrection(double jeteta, int runnumber){
// 	if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
// 		return 1.;
// 	}
// 	const double corrs[5] = {0.0, -0.454e-6, -0.952e-6, 1.378e-6, 0.0};
// 	const int run0 = 201000;
// 	double eta = fabs(jeteta);
// 	double corr = 0.;
// 	if(eta<1.3) corr = corrs[0];
// 	else if(eta<2.0) corr = corrs[1];
// 	else if(eta<2.5) corr = corrs[2];
// 	else if(eta<3.0) corr = corrs[3];
// 	else if(eta<5.0) corr = corrs[4];
// 	return (1.+corr*(runnumber-run0));
// }

// double Ntuple_Controller::JERCorrection(TLorentzVector jet, double dr, TString corr){
// 	double sf = jet.Pt();
// 	if (corr == "default") corr = jetCorrection;
// 	if(isData() || GetStrippedMCID()==DataMCType::DY_emu_embedded || GetStrippedMCID()==DataMCType::DY_mutau_embedded
// 			|| jet.Pt()<=10
// 			|| PFJet_matchGenJet(jet,dr)==TLorentzVector(0.,0.,0.,0.)
// 			){
// 		return sf;
// 	}else{
// 		double c = JetEnergyResolutionCorr(jet.Eta());
// 		if(corr.Contains("up")) c += JetEnergyResolutionCorrErr(jet.Eta());
// 		if(corr.Contains("down")) c -= JetEnergyResolutionCorrErr(jet.Eta());
// 		sf = std::max(0.,c*jet.Pt()+(1.-c)*PFJet_matchGenJet(jet,dr).Pt());
// 	}
// 	return sf;
// }

// TLorentzVector Ntuple_Controller::PFJet_matchGenJet(TLorentzVector jet, double dr){
// 	TLorentzVector genjet(0.,0.,0.,0.);
// 	for(unsigned int i=0;i<PFJet_NGenJetsNoNu();i++){
// 		if(PFJet_GenJetNoNu_p4(i).Vect()!=TVector3(0.,0.,0.)
// 				&& jet.DeltaR(PFJet_GenJetNoNu_p4(i))<dr
// 				){
// 			genjet = PFJet_GenJetNoNu_p4(i);
// 		}
// 	}
// 	return genjet;
// }

// double Ntuple_Controller::JetEnergyResolutionCorr(double jeteta){
// 	double eta = fabs(jeteta);
// 	double corr = 1.;
// 	if(eta<0.5) corr = 1.079;
// 	else if(eta<1.1) corr = 1.099;
// 	else if(eta<1.7) corr = 1.121;
// 	else if(eta<2.3) corr = 1.208;
// 	else if(eta<2.8) corr = 1.254;
// 	else if(eta<3.2) corr = 1.395;
// 	else if(eta<5.0) corr = 1.056;
// 	return corr;
// }

// double Ntuple_Controller::JetEnergyResolutionCorrErr(double jeteta){
// 	double eta = fabs(jeteta);
// 	double err = 0.;
// 	if(eta<0.5) err = 0.026;
// 	else if(eta<1.1) err = 0.027;
// 	else if(eta<1.7) err = 0.029;
// 	else if(eta<2.3) err = 0.046;
// 	else if(eta<2.8) err = 0.062;
// 	else if(eta<3.2) err = 0.063;
// 	else if(eta<5.0) err = 0.191;
// 	return err;
// }

/////////////////////////////////////////////////////////////////////
//
// Get jet four-vector
//
// Options:
//  - "run": corrects the jet pt to account for calorimeter degredation during data taking (only data).
//  - "JER": smears the jet pt in MC to match the resolution in data. Additional use of "up" or "down"
//           varies the correction by its uncertainty -> systematics
//  - "JEC": use this to estimate the impact of scale correction uncertainties on your result.
//           use "plus" for an upward variation. if you use nothing, the variation will be downward.
//

// TLorentzVector Ntuple_Controller::PFJet_p4(unsigned int i, TString corr){
// 	TLorentzVector vec = TLorentzVector(Ntp->PFJet_p4->at(i).at(1),Ntp->PFJet_p4->at(i).at(2),Ntp->PFJet_p4->at(i).at(3),Ntp->PFJet_p4->at(i).at(0));
// 	if (corr == "default") corr = jetCorrection;
// 	// apply run-dependent pT corrections
// 	if (corr.Contains("run")){
// 		vec.SetPerp(vec.Pt() * rundependentJetPtCorrection(vec.Eta(), RunNumber()));
// 	}
// 	if(corr.Contains("JER")){
// 		vec.SetPerp(JERCorrection(vec,0.25,corr));
// 	}
// 	if(corr.Contains("JEC")){
// 		if(corr.Contains("plus")) vec.SetPerp(vec.Pt() * (1 + PFJet_JECuncertainty(i)));
// 		else vec.SetPerp(vec.Pt() * (1 - PFJet_JECuncertainty(i)));
// 	}
// 	return vec;
// }

double Ntuple_Controller::TauSpinerGet(int SpinType){
#ifdef USE_TauSpinner
  if(!isData()){
    TauDecay taudecay;
    std::vector<SimpleParticle> tau_daughters, tau_daughters2;
    SimpleParticle tau, tau2;
    for(unsigned int i=0; i<NMCSignalParticles();i++){
      // check for signal Boson
      if(MCSignalParticle_pdgid(i)==25 || 
	 MCSignalParticle_pdgid(i)==36 || 
	 MCSignalParticle_pdgid(i)==22 || 
	 MCSignalParticle_pdgid(i)==23){
	if(MCSignalParticle_Tauidx(i).size()==2){
	  SimpleParticle X(MCSignalParticle_p4(i).Px(),MCSignalParticle_p4(i).Py(),MCSignalParticle_p4(i).Pz(),
			   MCSignalParticle_p4(i).E(),MCSignalParticle_pdgid(i));
	  bool tau1good(false),tau2good(false);
	  //first tau
	  unsigned int tauidx=MCSignalParticle_Tauidx(i).at(0);
	  Logger(Logger::Verbose) << "tau 1 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      int mypdgid=abs((int)MCTauandProd_pdgid(tauidx,t));
	      if(abs( mypdgid)==abs(PDGInfo::tau_plus) ){
		tau=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),MCTauandProd_p4(tauidx,t).Py(),MCTauandProd_p4(tauidx,t).Pz(),
				   MCTauandProd_p4(tauidx,t).E(),MCTauandProd_pdgid(tauidx,t));
	      }
	      else if((taudecay.isTauFinalStateParticle(mypdgid) && mypdgid!=PDGInfo::gamma)){
		tau_daughters.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),MCTauandProd_p4(tauidx,t).Py(),
						       MCTauandProd_p4(tauidx,t).Pz(),MCTauandProd_p4(tauidx,t).E(),
						       MCTauandProd_pdgid(tauidx,t)));
		tau1good=true;
	      }
	    }
	  }
	  // second tau
	  tauidx=MCSignalParticle_Tauidx(i).at(1);
	  Logger(Logger::Verbose) << "tau 2 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      int mypdgid=abs((int)MCTauandProd_pdgid(tauidx,t));
	      if(abs( mypdgid)==abs(PDGInfo::tau_plus)){
		    tau2=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
					MCTauandProd_p4(tauidx,t).Py(),
					MCTauandProd_p4(tauidx,t).Pz(),
					MCTauandProd_p4(tauidx,t).E(),
					MCTauandProd_pdgid(tauidx,t));
		    tau2good=true;
	      }
	      if((taudecay.isTauFinalStateParticle(mypdgid) && mypdgid!=PDGInfo::gamma)){
	    Logger(Logger::Verbose) << "isDaughter" << std::endl;
		if(tau_daughters.size()>0)tau2good=true;
		tau_daughters2.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
							MCTauandProd_p4(tauidx,t).Py(),
							MCTauandProd_p4(tauidx,t).Pz(),
							MCTauandProd_p4(tauidx,t).E(),
							MCTauandProd_pdgid(tauidx,t)));
	      }
	    }
	  }
	  if(tau1good && tau2good){
		Logger(Logger::Verbose)  << "Two Taus found: " << tau_daughters.size() << " " << tau_daughters2.size() << std::endl;
	    return TauSpinerInt.Get(SpinType,X,tau,tau_daughters,tau2,tau_daughters2);
	  }
	}
      }
    }
  }
#endif
  return 1.0;
}

/////////////////////////////////////////////////////////////////////
//
// Get tau four-vector
//
// Options:
//  - "scalecorr": corrects the tau energy scale depending on the decay mode (only MC and embedding).
//

// TLorentzVector Ntuple_Controller::PFTau_p4(unsigned int i, TString corr){
// 	TLorentzVector vec = TLorentzVector(Ntp->PFTau_p4->at(i).at(1),Ntp->PFTau_p4->at(i).at(2),Ntp->PFTau_p4->at(i).at(3),Ntp->PFTau_p4->at(i).at(0));
// 	if (corr == "default") corr = tauCorrection;
// 	if(!isData() || GetStrippedMCID() == DataMCType::DY_mutau_embedded){
// 		if(corr.Contains("scalecorr")){
// 			if(PFTau_hpsDecayMode(i)>0 && PFTau_hpsDecayMode(i)<5){
// 				vec *= 1.025+0.001*min(max(vec.Pt()-45.,0.),10.);
// 			}
// 			else if(PFTau_hpsDecayMode(i)>=10){
// 				vec *= 1.012+0.001*min(max(vec.Pt()-32.,0.),18.);
// 			}
// 		}
// 		if(corr.Contains("met")){
// 			if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.03);
// 			else vec.SetPerp(vec.Perp() * 0.97);
// 		}
// 	}
// 	return vec;
// }


bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &tau_idx){
  for(unsigned int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==(int)parent_pdgid){
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	if((unsigned int)MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
	  Logger(Logger::Warning) << "INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
	  return false;
	}
      }
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	unsigned int tauidx=MCSignalParticle_Tauidx(i).at(j);
	Logger(Logger::Verbose) << "MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << " " << Ntp->MCTau_JAK->size() << " " << Ntp->MCTauandProd_pdgid->size() << std::endl;
	if((int)MCTau_JAK(tauidx)==tau_jak){ tau_idx=tauidx;Boson_idx=i;return true;}
      }
    }
  }
  return false;
}

bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx){
  for(unsigned int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==parent_pdgid){
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
        if((unsigned int)MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
          Logger(Logger::Warning) << "INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
          return false;
        }
      }
      if(MCSignalParticle_Tauidx(i).size()==2){
	tau1_idx=MCSignalParticle_Tauidx(i).at(0);
        tau2_idx=MCSignalParticle_Tauidx(i).at(1);
	Boson_idx=i;
	return true;
      }
    }
  }
  return false;
}





  int Ntuple_Controller::getBitOfGivenTrigger(TString tname){

    TFile *currentFile = Ntp->fChain->GetCurrentFile();
    TH1F  *hLLRCounters = (TH1F*)currentFile->Get("HTauTauTree/Counters");

    int ibit=-1;
    for(unsigned int iBinX=4;iBinX<(unsigned int)(hLLRCounters->GetNbinsX()+1);++iBinX){
      TString name = hLLRCounters->GetXaxis()->GetBinLabel(iBinX);
      if(name.Contains(tname))ibit =  iBinX-4;
    }
    return ibit;
  }
  

bool Ntuple_Controller::GetTriggerIndex(TString n,  int &i){
  for(i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      if(name.Contains(n))return true;
    } 
	return false;
 }


bool  Ntuple_Controller::CheckIfAnyPassed(  std::vector<int> list){
  for(unsigned int itrig = 0; itrig < list.size(); itrig++){
    if(TriggerAccept(list.at(itrig)))      return true;
    
  }
  return false;
}


std::vector<int> Ntuple_Controller::GetVectorTriggers(TString n){
    std::vector<int> out;
    for(int i=0; i<NTriggers();i++){
	TString name=TriggerName(i);
	if(name.Contains(n)) out.push_back(i) ;
      } 
	  return out;
 }

std::vector<int> Ntuple_Controller::GetVectorTriggers(std::vector<TString>  v){
    std::vector<int> out;
    for(int i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      for(unsigned int j=0; j<v.size(); j++){
	if(name.Contains(v.at(j))) out.push_back(i) ;
      }
    } 
    return out;
}

std::vector<int> Ntuple_Controller::GetVectorTriggersFullMatch(std::vector<TString>  v){
    std::vector<int> out;
    for(int i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      bool cpattern(true);
      for(unsigned int j=0; j<v.size(); j++){
	if(!name.Contains(v.at(j))) cpattern =false;
      }
      if(cpattern)  out.push_back(i) ;
    } 
    return out;
}

std::vector<int> Ntuple_Controller::GetVectorCrossTriggers(TString n1,TString n2,TString f1,TString f2){
    std::vector<int> out;



    for(int i=0; i<NTriggers();i++){
      TString name=TriggerName(i);
      if(name.Contains(n1) &&  name.Contains(n2)  && (!name.Contains(f1)  && !name.Contains(f2))   ) out.push_back(i) ;
    } 
    return out;
}



// // PFTau significance, using the reffited primary and secondary vertices
// double Ntuple_Controller::PFTau_FlightLength_significance(unsigned int i) {
// 	TVector3 PVpos = PFTau_TIP_primaryVertex_pos(i);
// 	TMatrixTSym<double> PVcov = PFTau_TIP_primaryVertex_cov(i);
// 	TVector3 SVpos = PFTau_TIP_secondaryVertex_pos(i);
// 	TMatrixTSym<double> SVcov = PFTau_TIP_secondaryVertex_cov(i);

// 	return PFTau_FlightLength_significance(PVpos, PVcov, SVpos, SVcov);
// }

// // calculate flight length significance from primary and secondary vertex info
 double Ntuple_Controller::PFTau_FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
   TVector3 SVPV = sv - pv;
   TVectorF FD;
   FD.ResizeTo(3);
   FD(0) = SVPV.X();
   FD(1) = SVPV.Y();
   FD(2) = SVPV.Z();

   TMatrixT<double> PVcv;
   PVcv.ResizeTo(3,3);
   for(int nr =0; nr<PVcov.GetNrows(); nr++){
     for(int nc =0; nc<PVcov.GetNcols(); nc++){
       PVcv(nr,nc) = PVcov(nr,nc);
     }
   }
   TMatrixT<double> SVcv;
   SVcv.ResizeTo(3,3);
   for(int nr =0; nr<SVcov.GetNrows(); nr++){
     for(int nc =0; nc<SVcov.GetNcols(); nc++){
       SVcv(nr,nc) = SVcov(nr,nc);
     }
   }

   TMatrixT<double> SVPVMatrix(3,1);
   for(int i=0; i<SVPVMatrix.GetNrows();i++){
     SVPVMatrix(i,0)=FD(i);
   }

   TMatrixT<double> SVPVMatrixT=SVPVMatrix;
   SVPVMatrixT.T();

   TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
   double sigmaabs = sqrt(lambda2(0,0))/SVPV.Mag();
   double sign = SVPV.Mag()/sigmaabs;

   return sign;
 }

//// Generator Information
int Ntuple_Controller::matchTruth(TLorentzVector tvector){
	double testdr=0.3;
	int pdgid = 0;
	for(unsigned int i=0;i<NMCParticles();i++){
		if(MCParticle_p4(i).Pt()>0.){
			if(tvector.DeltaR(MCParticle_p4(i))<testdr){
				testdr = tvector.DeltaR(MCParticle_p4(i));
				pdgid = MCParticle_pdgid(i);
			}
		}
	}
	return pdgid;
}
bool Ntuple_Controller::matchTruth(TLorentzVector tvector, int pid, double dr){
	if (getMatchTruthIndex(tvector, pid, dr) >= 0) return true;
	return false;
}
int Ntuple_Controller::getMatchTruthIndex(TLorentzVector tvector, int pid, double dr){
	int index = -9;
	for(unsigned int i=0;i<NMCParticles();i++){
		if(MCParticle_p4(i).Pt()>0.){
			if(fabs(MCParticle_pdgid(i))==pid){
				if(tvector.DeltaR(MCParticle_p4(i))<dr) index = i;
			}
		}
	}
	return index;
}

int Ntuple_Controller::MCTau_true3prongAmbiguity(unsigned int i){
	int amb = -9;

	int j_a1 = MCTau_getDaughterOfType(i, PDGInfo::a_1_plus, true);
	if (j_a1 < 0){
		Logger(Logger::Warning) << "MCTau has no a1 decay product." << std::endl;
		return -9;
	}

	TLorentzVector genTauh	= MCTau_p4(i);
	TLorentzVector genA1	= MCTauandProd_p4(i, j_a1);

	TLorentzVector GenA1_boosted(BoostToRestFrame(genTauh, genA1));
	double dotProduct = GenA1_boosted.Vect().Dot( genTauh.Vect() );
	if(dotProduct < 0){
		amb = MultiProngTauSolver::plus;
	}
	else if(dotProduct > 0){
		amb = MultiProngTauSolver::minus;
	}
	else if(dotProduct == 0){
		amb = MultiProngTauSolver::zero;
	}

	return amb;
}
int Ntuple_Controller::MCTau_getDaughterOfType(unsigned int i_mcTau, int daughter_pdgid, bool ignoreCharge /*= true*/){
	int matchedIndex = -1;
	for(int i_dau=0; i_dau < NMCTauDecayProducts(i_mcTau); i_dau++){
		if( ignoreCharge ){
			if( abs(MCTauandProd_pdgid(i_mcTau, i_dau)) == abs(daughter_pdgid) )
				matchedIndex = i_dau;
		}
		else{
			if( MCTauandProd_pdgid(i_mcTau, i_dau) == daughter_pdgid )
				matchedIndex = i_dau;
		}
	}
	return matchedIndex;
}
// int Ntuple_Controller::matchTauTruth(unsigned int i_hpsTau, bool onlyHadrDecays /*= false*/){
// 	int matchedIndex = -1;
// 	double minDr = 999;
// 	for(int i=0; i < NMCTaus(); i++){
// 		if( onlyHadrDecays && (MCTau_JAK(i) <= 2) ) continue; // exclude decays to electrons and muons
// 		double dr = MCTau_p4(i).DeltaR(PFTau_p4(i_hpsTau));
// 		if( dr < 0.5 && dr < minDr ){
// 			matchedIndex = i;
// 			minDr = dr;
// 		}
// 	}
// 	return matchedIndex;
// }

//gets two TLorentzVectors (1. defines RF/Boost, 2. is boosted), creates a copy of the second and boosts it
TLorentzVector Ntuple_Controller::BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2){
	TVector3 boostvector = TLV1.BoostVector();
	TLorentzVector boosted_TLV2(TLV2);
	boosted_TLV2.Boost(-boostvector);
	return boosted_TLV2;
}

// get visible/invisible part of gen Tau 4-vector
TLorentzVector Ntuple_Controller::MCTau_invisiblePart(unsigned int i){
	TLorentzVector lv(0,0,0,0);
	for(int i_dau=1; i_dau < NMCTauDecayProducts(i); i_dau++){
		if( abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_e ||
			abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_mu ||
			abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_tau)
				lv += MCTauandProd_p4(i, i_dau);
	}
	if (lv.Pt() == 0.0)
		Logger(Logger::Warning) << "This tau decay has no neutrinos!?" << std::endl;

	return lv;
}
TLorentzVector Ntuple_Controller::MCTau_visiblePart(unsigned int i){
	TLorentzVector lv = MCTau_p4(i);
	lv -= MCTau_invisiblePart(i);
	return lv;
}

//// draw decay chain
void Ntuple_Controller::printMCDecayChainOfEvent(bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	Logger(Logger::Info) << "=== Draw MC decay chain of event ===" << std::endl;
	for(unsigned int par = 0; par < NMCParticles(); par++){
		if ( !MCParticle_hasMother(par) )
			printMCDecayChainOfMother(par, printStatus, printPt, printEtaPhi, printQCD);
	}
}
void Ntuple_Controller::printMCDecayChainOfMother(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	Logger(Logger::Info) << "Draw decay chain for mother particle at index " << par << " :" << std::endl;
	printMCDecayChain(par, 0, printStatus, printPt, printEtaPhi, printQCD);
}
void Ntuple_Controller::printMCDecayChain(unsigned int par, unsigned int level, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	std::ostringstream out;
	for(unsigned int l = 0; l < level; l++) out << "    ";
	out << MCParticleToString(par, printStatus, printPt, printEtaPhi);
	if (!printQCD && MCParticle_pdgid(par) == 92) {out << "    QCD stuff truncated ...\n"; std::cout << out.str(); return;}
	out << "\n";
	std::cout << out.str();
	for (unsigned int i_dau = 0; i_dau < MCParticle_childidx(par).size(); i_dau++){
		if (MCParticle_childidx(par).at(i_dau) != (int)par) // skip cases where particles are daughters of themselves
			printMCDecayChain(MCParticle_childidx(par).at(i_dau), level+1, printStatus, printPt, printEtaPhi, printQCD);
	}
}

std::string Ntuple_Controller::MCParticleToString(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi){
	std::ostringstream out;
	out << "+-> ";
	out << PDGInfo::pdgIdToName( MCParticle_pdgid(par) );
	if (printStatus) out << " (status = " << MCParticle_status(par) << ", idx = " << par << ")";
	if (printPt || printEtaPhi) out <<  " [";
	if (printPt) out << "pT = " << MCParticle_p4(par).Pt() << "GeV";
	if (printEtaPhi) out << " eta = " << MCParticle_p4(par).Eta() << " phi = " << MCParticle_p4(par).Phi();
	if (printPt || printEtaPhi) out << "]";
	return out.str();
}
bool Ntuple_Controller::CheckDecayID(unsigned  int jak1, unsigned int jak2){
  // if(id!=998) return false;
  bool decayid= false;
  for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
    // std::cout<<"NC: Ntp->MCSignalParticle_p4->size()  "<<Ntp->MCSignalParticle_p4->size()<<std::endl;
    // std::cout<<"NC: Ntp->MCSignalParticle_Tauidx->at(iz).size()  "<<Ntp->MCSignalParticle_Tauidx->at(iz).size()<<std::endl;
    // std::cout<<"  MCSignalParticle_pdgid->at(iz)  "<< Ntp-> MCSignalParticle_pdgid->at(iz) <<std::endl;
    if(fabs(Ntp-> MCSignalParticle_pdgid->at(iz) )!=23) return false;
    if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
      //     std::cout<<" Ntp->MCTau_JAK-> "<< Ntp->MCTau_JAK->size()<<std::endl;
      if(Ntp->MCTau_JAK->at(0) == jak1 and Ntp->MCTau_JAK->at(1) ==jak2 ){ decayid = true;}
      else if(Ntp->MCTau_JAK->at(0) ==jak2  and Ntp->MCTau_JAK->at(1) ==jak1){decayid  = true;}
    }
  }
  return decayid;
}
TLorentzVector Ntuple_Controller::GetTruthTauLV(unsigned int jak){
  TLorentzVector tau(0,0,0,0);
  bool DecayOK = false;
  unsigned int tauIndex;
  for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
    if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
      if(Ntp->MCTau_JAK->at(0) == jak){tauIndex=0; DecayOK = true;}
      else if(  Ntp->MCTau_JAK->at(1) ==jak ){ tauIndex=1; DecayOK = true;}
      if(DecayOK){
        tau = TLorentzVector(Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(1),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(2),
                             Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(3),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(0));

      }
    }
  }
  return tau;
}


TLorentzVector Ntuple_Controller::GetTruthTauProductLV(unsigned int jak, int pdgID){
  TLorentzVector tauProd(0,0,0,0);
  bool DecayOK = false;
  unsigned int tauIndex;
  for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
    if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
      if(Ntp->MCTau_JAK->at(0) == jak){tauIndex=0; DecayOK = true;}
      else if(  Ntp->MCTau_JAK->at(1) ==jak ){ tauIndex=1; DecayOK = true;}
      if(DecayOK){
        unsigned int NDec;
        if(0<=tauIndex && tauIndex<NMCTaus()){ NDec = Ntp->MCTauandProd_p4->at(tauIndex).size();}
        else NDec= 0;
        for(unsigned int iProd =0; iProd < NDec; iProd++ ){
          if(abs( Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd))==pdgID){
            //if( Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd)==pdgID){
              tauProd = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(2),
                                       Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(0));
              // }
          }
        }
      }
    }
  }
    return tauProd;
}


std::vector<TLorentzVector> Ntuple_Controller::GetTruthPionsFromA1(){
        TLorentzVector SSPion1(0,0,0,0);
        TLorentzVector SSPion2(0,0,0,0);
        TLorentzVector OSPion(0,0,0,0);
        unsigned int OSMCPionIndex;
        unsigned int SSMCPion1Index;
        unsigned int SSMCPion2Index;

        std::vector<TLorentzVector>  output;

        bool DecayOK = false;
        unsigned int tauIndex;
        for(unsigned int iz =0; iz<Ntp->MCSignalParticle_p4->size(); iz++){
          if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
            if(Ntp->MCTau_JAK->at(0) == 5){tauIndex=0; DecayOK = true;}
            else if(  Ntp->MCTau_JAK->at(1) ==5 ){ tauIndex=1; DecayOK = true;}
            if(DecayOK){
              unsigned int NDec;
              if(0<=tauIndex && tauIndex<NMCTaus()){ NDec = Ntp->MCTauandProd_p4->at(tauIndex).size();}
              else NDec= 0;
              int nplus =0, nminus = 0;
              for(unsigned int iProd =0; iProd < NDec; iProd++ ){
		// if(abs(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd))==211){	  std::cout<<" pdgId     "<< Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd)<< std::endl;
		//   std::cout<< " print mc again   "<<  Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(1) << "  "<<Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(2)	  << "  "<<Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(3)	  << "  "<<Ntp->MCTauandProd_p4->at(tauIndex).at(iProd).at(0)	<<	  std::endl; }

                if( Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd)== 211) nplus++;
                if( Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd)==-211) nminus++;
              }

              //    std::cout<<" nplus nminus  "<< nplus << "  "<< nminus<<std::endl;

              if(nplus == 1 && nminus==2){
                int nss=0;
                for(unsigned int iProd1 =0; iProd1 < NDec; iProd1++ ){
                  if(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd1)== 211){
                    OSPion = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(2), Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(0));
                    OSMCPionIndex=iProd1;
                  }

                  if(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd1)==-211 && nss ==0){
                    nss++;
                    SSPion1 = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(2), Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(0));
                    SSMCPion1Index=iProd1;
                  }
                  //              std::cout<<" nss "<< nss << " iProd1 "<<iProd1<<std::endl;
                  if(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd1)==-211 &&  nss == 1){
                    SSPion2 = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(2), Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd1).at(0));
                    SSMCPion2Index=iProd1;
                  }
                }
              }

              if(nplus == 2 && nminus==1){
                int nss=0;
                for(unsigned int iProd2 =0; iProd2 < NDec; iProd2++ ){
                  if(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd2)== -211){
                    OSPion = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(2), Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(0));
                    OSMCPionIndex=iProd2;
                  }

                  if(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd2)==211 && nss ==0){
                    nss++;
                    SSPion1 = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(2), Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(0));
                    SSMCPion1Index=iProd2;
                  }
                  if(Ntp->MCTauandProd_pdgid->at(tauIndex).at(iProd2)==211 && nss ==1){
                    SSPion2 = TLorentzVector(Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(1),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(2), Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(3),Ntp->MCTauandProd_p4->at(tauIndex).at(iProd2).at(0));
                    SSMCPion2Index=iProd2;
                  }
                }
              }
            }
          }
        }

        output.push_back(OSPion);
        output.push_back(SSPion1);
        output.push_back(SSPion2);
        return output;
}





//// Trigger Information
// bool Ntuple_Controller::TriggerAccept(TString n){
//   unsigned int i=0;
//   if(GetTriggerIndex(n,i))return TriggerAccept(i);
//   return false;
// }

// unsigned int Ntuple_Controller::HLTPrescale(TString n){
//   unsigned int i=0;
//   if(GetTriggerIndex(n,i))return HLTPrescale(i);
//   return 1;
// }

// unsigned int Ntuple_Controller::L1SEEDPrescale(TString n){
//   unsigned int i=0;
//   if(GetTriggerIndex(n,i))return L1SEEDPrescale(i);
//   return 1;
// }




// double Ntuple_Controller::matchTrigger(TLorentzVector obj, std::vector<TString> trigger, std::string objectType){
//   unsigned int id = 0;
//   TLorentzVector triggerObj(0.,0.,0.,0.);
//   if(objectType=="tau"){
//     id = 84;
//   }
//   if(objectType=="muon"){
//     id = 83;
//   }
//   if(objectType=="electron"){
//     id = 82;
//   }
  
//   double minDR = 100.;
//   for(unsigned int i_trig = 0; i_trig < trigger.size(); i_trig++){
//     for(int i=0;i<NHLTTrigger_objs();i++){
//       if(HLTTrigger_objs_trigger(i).find(trigger.at(i_trig)) != string::npos){
// 	for(int j=0;j<NHLTTrigger_objs(i);j++){
// 	  if(HLTTrigger_objs_Id(i,j)==(int)id){
// 	    triggerObj.SetPtEtaPhiE(HLTTrigger_objs_Pt(i,j),
// 				    HLTTrigger_objs_Eta(i,j),
// 				    HLTTrigger_objs_Phi(i,j),
// 				    HLTTrigger_objs_E(i,j));
// 	  }
// 	  if( triggerObj.Pt()>0. && obj.Pt()>0. ) {
// 	    double dr = obj.DeltaR(triggerObj);
// 	    if (dr < minDR) minDR = dr;
// 	  }
// 	}
//       }
//     }
//   }
//   return minDR;
// }
// bool Ntuple_Controller::matchTrigger(TLorentzVector obj, double dr_cut, std::vector<TString> trigger, std::string objectType){
// 	double dr = matchTrigger(obj, trigger, objectType);
// 	return dr < dr_cut;
// }
// bool Ntuple_Controller::matchTrigger(TLorentzVector obj, double dr_cut, TString trigger, std::string objectType){
// 	std::vector<TString> triggerVec;
// 	triggerVec.push_back(trigger);
// 	return matchTrigger(obj, dr_cut, triggerVec, objectType);
// }

bool Ntuple_Controller::isPVCovAvailable(){ // sometimes returns zero size matrix (rare)
  if(Ntp->pv_cov->size()!=6)  return false; 
  return true;

}
TMatrixTSym<float> Ntuple_Controller::PFTau_TIP_primaryVertex_cov(){
  TMatrixTSym<float> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(int k=j;k<LorentzVectorParticle::NVertex;k++){
      //if(j==k) V_cov(i,j)=pow(0.0001,2.0);
      V_cov(j,k)=Ntp->pv_cov->at(l);
      V_cov(k,j)=Ntp->pv_cov->at(l);
      l++;
    }
  }
  //  std::cout<<"  pv is good"<< std::endl; V_cov.Print();
  return  V_cov;
}

TMatrixTSym<double> Ntuple_Controller::PFTau_TIP_secondaryVertex_cov(unsigned int i){
  TMatrixTSym<double> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(int k=j;k<LorentzVectorParticle::NVertex;k++){
      V_cov(j,k)=Ntp->PFTauSVCov->at(i).at(l);
      V_cov(k,j)=Ntp->PFTauSVCov->at(i).at(l);
      l++;
    }
  }
  //    std::cout<<"  sv is good"<< std::endl; V_cov.Print();
  return  V_cov;
}

 LorentzVectorParticle Ntuple_Controller::PFTau_a1_lvp(unsigned int i){
   TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
   TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
   int l=0;
   if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
     for(int k=0; k<LorentzVectorParticle::NLorentzandVertexPar; k++){
       a1_par(k,0)=Ntp->PFTau_a1_lvp->at(i).at(k);
       for(int j=k; j<LorentzVectorParticle::NLorentzandVertexPar; j++){
 	a1_cov(k,j)=Ntp->PFTau_a1_cov->at(i).at(l);
 	a1_cov(j,k)=Ntp->PFTau_a1_cov->at(i).at(l);
 	l++;
       } 
     }
   }
   return LorentzVectorParticle(a1_par,a1_cov,Ntp->PFTau_a1_pdgid->at(i),Ntp->PFTau_a1_charge->at(i),Ntp->PFTau_a1_B->at(i));
 }

// std::vector<TrackParticle> Ntuple_Controller::PFTau_daughterTracks(unsigned int i){
//   std::vector<TrackParticle> daughter;
//   for(unsigned int d=0;d<Ntp->PFTau_daughterTracks_poca->at(i).size();d++){
//     TMatrixT<double>    a1_par(TrackParticle::NHelixPar,1);
//     TMatrixTSym<double> a1_cov(TrackParticle::NHelixPar);
//     int l=0;
//     for(int k=0; k<TrackParticle::NHelixPar; k++){
//       a1_par(k,0)=Ntp->PFTau_daughterTracks->at(i).at(d).at(k);
//       for(int j=k; j<TrackParticle::NHelixPar; j++){
// 	a1_cov(k,j)=Ntp->PFTau_daughterTracks->at(i).at(d).at(l);
// 	a1_cov(j,k)=Ntp->PFTau_daughterTracks->at(i).at(d).at(l);
// 	l++;
//       }
//     }
//     daughter.push_back(TrackParticle(a1_par,a1_cov,Ntp->PFTau_daughterTracks_pdgid->at(i).at(d),Ntp->PFTau_daughterTracks_M->at(i).at(d),Ntp->PFTau_daughterTracks_charge->at(i).at(d),Ntp->PFTau_daughterTracks_B->at(i).at(d)));
//   }
//   return daughter;
// }

// std::vector<TVector3> Ntuple_Controller::PFTau_daughterTracks_poca(unsigned int i){
//   std::vector<TVector3> poca;
//   for(unsigned int k=0;k<Ntp->PFTau_daughterTracks_poca->at(i).size();k++){
//     poca.push_back(TVector3(Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(0),Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(1),Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(2)));
//   }
//   return poca;
// }
   
double Ntuple_Controller::dxySigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return (-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt();
}
double Ntuple_Controller::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(dxySigned(fourvector, poca, vtx));
}


double Ntuple_Controller::dzSigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2);
}
double Ntuple_Controller::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(dzSigned(fourvector, poca, vtx));
}

// double Ntuple_Controller::vertexSignificance(TVector3 vec, unsigned int vertex){
// 	if(vertex>=0 && vertex<NVtx()){
// 		const double elm[3] = {(vec.X()-Vtx(vertex).X()),(vec.Y()-Vtx(vertex).Y()),(vec.Z()-Vtx(vertex).Z())};
// 		TVectorD diff(3,elm);
// 		TMatrixF M(Vtx_Cov(vertex));
// 		if(M.IsValid()){
// 			double mag = diff.Norm2Sqr();
// 			double sim = M.Similarity(diff);
// 			return mag/sqrt(sim);
// 		}
// 	}
// 	return 999;
// }

// check if given lepton was used for MVA-MET calculation
// bool Ntuple_Controller::findCorrMVASrcMuon(unsigned int muon_idx, int &mvaSrcMuon_idx, float &dR ){
// 	float minDr = 1000;
// 	float dr = 1001;
// 	for(unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVA_srcMuons(); i_mvaLep++){
// 		 dr = Tools::dr(Muon_p4(muon_idx),MET_CorrMVA_srcMuon_p4(i_mvaLep));
// 		 if ((dr < 0.05) && (dr < minDr)) {
// 			 minDr = dr;
// 			 dR = dr;
// 			 mvaSrcMuon_idx = i_mvaLep;
// 		 }
// 	}
// 	if (minDr < 0.05) return true;
// 	else return false;
// }
// bool Ntuple_Controller::findCorrMVASrcElectron(unsigned int elec_idx, int &mvaSrcElectron_idx, float &dR ){
// 	float minDr = 1000;
// 	float dr = 1001;
// 	for (unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVA_srcElectrons(); i_mvaLep++){
// 		 dr = Tools::dr(Electron_p4(elec_idx),MET_CorrMVA_srcElectron_p4(i_mvaLep));
// 		 if ((dr < 0.05) && (dr < minDr)) {
// 			 minDr = dr;
// 			 dR = dr;
// 			 mvaSrcElectron_idx = i_mvaLep;
// 		 }
// 	}
// 	if (minDr < 0.05) return true;
// 	else return false;
// }
// bool Ntuple_Controller::findCorrMVASrcTau(unsigned int tau_idx, int &mvaSrcTau_idx, float &dR ){
// 	float minDr = 1000;
// 	float dr = 1001;
// 	for(unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVA_srcTaus(); i_mvaLep++){
// 		 dr = Tools::dr(PFTau_p4(tau_idx),MET_CorrMVA_srcTau_p4(i_mvaLep));
// 		 if ((dr < 0.05) && (dr < minDr)) {
// 			 minDr = dr;
// 			 dR = dr;
// 			 mvaSrcTau_idx = i_mvaLep;
// 		 }
// 	}
// 	if (minDr < 0.05) return true;
// 	else return false;
// }
// bool Ntuple_Controller::findCorrMVAMuTauSrcMuon(unsigned int muon_idx, int &mvaMuTauSrcMuon_idx, float &dR ){
// 	float minDr = 1000;
// 	float dr = 1001;
// 	for (unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVAMuTau_srcMuons(); i_mvaLep++){
// 		 dr = Tools::dr(Muon_p4(muon_idx),MET_CorrMVAMuTau_srcMuon_p4(i_mvaLep));
// 		 if ((dr < 0.05) && (dr < minDr)) {
// 			 minDr = dr;
// 			 dR = dr;
// 			 mvaMuTauSrcMuon_idx = i_mvaLep;
// 		 }
// 	}
// 	if (minDr < 0.05) return true;
// 	else return false;
// }
// bool Ntuple_Controller::findCorrMVAMuTauSrcTau(unsigned int tau_idx, int &mvaMuTauSrcTau_idx, float &dR ){
// 	float minDr = 1000;
// 	float dr = 1001;
// 	for (unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVAMuTau_srcTaus(); i_mvaLep++){
// 		 dr = Tools::dr(PFTau_p4(tau_idx),MET_CorrMVAMuTau_srcTau_p4(i_mvaLep));
// 		 if ((dr < 0.05) && (dr < minDr)) {
// 			 minDr = dr;
// 			 dR = dr;
// 			 mvaMuTauSrcTau_idx = i_mvaLep;
// 		 }
// 	}
// 	if (minDr < 0.05) return true;
// 	else return false;
// }

// function to sort any objects by any value in descending order
std::vector<int> Ntuple_Controller::sortObjects(std::vector<int> indices, std::vector<double> values){
	if (indices.size() != values.size()){
		Logger(Logger::Warning) << "Please make sure indices and values have same size for sorting. Abort." << std::endl;
		return std::vector<int>();
	}
	// create vector of pairs to allow for sorting by value
	std::vector< std::pair<int, double> > pairs;
	for(unsigned int i = 0; i<values.size(); i++ ){
		pairs.push_back( std::make_pair(indices.at(i),values.at(i)) );
	}
	// sort vector of pairs
	std::sort(pairs.begin(), pairs.end(), sortIdxByValue());
	// create vector of indices in correct order
	std::vector<int> sortedIndices;
	for(unsigned int i = 0; i<pairs.size(); i++){
		sortedIndices.push_back(pairs.at(i).first);
	}
	return sortedIndices;
}

// std::vector<int> Ntuple_Controller::sortDefaultObjectsByPt(TString objectType){
// 	std::vector<int> indices;
// 	std::vector<double> values;
// 	if (objectType == "Jets" || objectType == "PFJets"){
// 	  for(unsigned int i = 0; i<NPFJets(); i++ ){
// 	    indices.push_back(i);
// 	    values.push_back(PFJet_p4(i).Pt());
// 	  }
// 	}
// 	else if(objectType == "Taus" || objectType == "PFTaus"){
// 	  for(unsigned int i = 0; i<NPFTaus(); i++ ){
// 	    indices.push_back(i);
// 	    values.push_back(PFTau_p4(i).Pt());
// 	  }
// 	}
// 	else if(objectType == "Muons"){
// 	  for (unsigned int i = 0; i<NMuons(); i++ ){
// 	    indices.push_back(i);
// 	    values.push_back(Muon_p4(i).Pt());
// 	  }
// 	}
// 	else if(objectType == "Electrons"){
// 	  for (unsigned int i = 0; i<NElectrons(); i++ ){
// 	    indices.push_back(i);
// 	    values.push_back(Electron_p4(i).Pt());
// 	  }
// 	}
// 	else{
// 	  Logger(Logger::Warning) << "sortDefaultObjectsByPt is only implemented for Jets, Taus, Muons and Electrons. Abort." << std::endl;
// 	  return std::vector<int>();
// 	}
// 	return sortObjects(indices, values);
// }



// obtain, or create and store, SVFit results from/on dCache
#ifdef USE_SVfit
SVFitObject* Ntuple_Controller::getSVFitResult_MuTauh(SVFitStorage& svFitStor, TString metType, unsigned muIdx, unsigned tauIdx, unsigned rerunEvery /* = 5000 */, TString suffix /* ="" */, double scaleMu /* =1 */, double scaleTau /* =1 */) {
	 // configure svfitstorage on first call
	if ( !svFitStor.isConfigured() ) svFitStor.Configure(GetInputDatasetName(), suffix);
	// get SVFit result from cache
	SVFitObject* svfObj = svFitStor.GetEvent(RunNumber(), LuminosityBlock(), EventNumber());
	// if obtained object is not valid, create and store it
	if (!svfObj->isValid()) {
		runAndSaveSVFit_MuTauh(svfObj, svFitStor, metType, muIdx, tauIdx, scaleMu, scaleTau);
	}
	else{
		// calculate every N'th event and compare with what is stored
		if( (EventNumber() % rerunEvery) == 123){
			SVFitObject* newSvfObj = new SVFitObject();
			runAndSaveSVFit_MuTauh(newSvfObj, svFitStor, metType, muIdx, tauIdx, scaleMu, scaleTau, false); // will not be saved in output files

			if (*svfObj == *newSvfObj){
				Logger(Logger::Info) << "Recalculation of SVFit object gave same result." << std::endl;
			}
			else {
				Logger(Logger::Warning) << "Recalculation of SVFit object gave DIFFERENT result!!" <<
				"\n\told: mass = " << svfObj->get_mass() << " +/- " << svfObj->get_massUncert() << ", pt = " << svfObj->get_pt() << " +/- " << svfObj->get_ptUncert() <<
				"\n\tnew: mass = " << newSvfObj->get_mass() << " +/- " << newSvfObj->get_massUncert() << ", pt = " << newSvfObj->get_pt() << " +/- " << newSvfObj->get_ptUncert() <<
				"\n\tSmall discrepancies could be caused by fitting details. It's up to you whether to ignore them." << std::endl;
			}
			delete newSvfObj;
		}
	}
	return svfObj;
}
SVFitObject* Ntuple_Controller::getSVFitResult_MuTau3p(SVFitStorage& svFitStor, TString metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, TString suffix /* ="" */, double scaleMu /* =1 */, double scaleTau /* =1 */) {
	 // configure svfitstorage on first call
	if ( !svFitStor.isConfigured() ) svFitStor.Configure(GetInputDatasetName(), suffix);
	// get SVFit result from cache
	SVFitObject* svfObj = svFitStor.GetEvent(RunNumber(), LuminosityBlock(), EventNumber());
	// if obtained object is not valid, create and store it
	if (!svfObj->isValid()) {
		runAndSaveSVFit_MuTau3p(svfObj, svFitStor, metType, muIdx, tauLV, neutrino, scaleMu, scaleTau);
	}
	return svfObj;
}

// create SVFitObject from standard muon and standard tau_h
void Ntuple_Controller::runAndSaveSVFit_MuTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, unsigned tauIdx, double scaleMu, double scaleTau, bool save /*= true*/) {
	objects::MET met(this, metType);
	SVfitProvider svfProv(this, met, "Mu", muIdx, "Tau", tauIdx, 1, scaleMu, scaleTau);
	*svfObj = svfProv.runAndMakeObject();
	if (svfObj->isValid()) {
		// store only if object is valid
		if (save) svFitStor.SaveEvent(RunNumber(), LuminosityBlock(), EventNumber(), svfObj);
	} else {
		Logger(Logger::Error) << "Unable to create a valid SVFit object." << std::endl;
	}
}
// create SVFitObject from standard muon and fully reconstructed 3prong tau
void Ntuple_Controller::runAndSaveSVFit_MuTau3p(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, double scaleMu, double scaleTau, bool save /*= true*/){
	objects::MET met(this, metType);
	met.subtractNeutrino(neutrino);
	SVfitProvider svfProv(this, met, "Mu", muIdx, tauLV, 1, scaleMu, scaleTau);
	*svfObj = svfProv.runAndMakeObject();
	if (svfObj->isValid()) {
		// store only if object is valid
		if (save) svFitStor.SaveEvent(RunNumber(), LuminosityBlock(), EventNumber(), svfObj);
	} else {
		Logger(Logger::Error) << "Unable to create a valid SVFit object." << std::endl;
	}
}




#endif // USE_SVfit
