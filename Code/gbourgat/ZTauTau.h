#ifndef ZTauTau_h
#define ZTauTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
//#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "ReferenceScaleFactors.h"
#include "ScaleFactor.h"
#include "Objects.h"
#include "PUReweight.h"
#include "tauTrigSFreader.h"
#include "DataMCCorrections.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

class ZTauTau : public Selection {

 public:
  ZTauTau(TString Name_, TString id_);
  virtual ~ZTauTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {Trigger=0,
	     Id_and_Kin,
	     NPairsFound,
	     Tau1Isolation,
	     Tau2Isolation,
	     LeptonVeto,
	     PairCharge,
	     PairMass,
	     //MTM,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  DataMCCorrections DataMC_Corr;
  tauTrigSFreader tauTrgSF;

  ClassicSVfit svfitAlgo1;
  //ClassicSVfit svfitAlgo2;
  //  SVFitStorage svfitstorage;


 private:
  // Selection Variables and Histos
  
  std::vector<TH1D> Tau1PT;
  std::vector<TH1D> Tau1E;
  std::vector<TH1D> Tau1Mass;
  std::vector<TH1D> Tau1Phi;
  std::vector<TH1D> Tau1Eta;
  std::vector<TH1D> Tau1dz;

  std::vector<TH1D> Tau2PT;
  std::vector<TH1D> Tau2E;
  std::vector<TH1D> Tau2Mass;
  std::vector<TH1D> Tau2Phi;
  std::vector<TH1D> Tau2Eta;
  std::vector<TH1D> Tau2dz;
  /*
  std::vector<TH1D> Tau1isolation;
  std::vector<TH1D> Tau2isolation;
  
  std::vector<TH1D> againstElectronVLooseMVA6_Tau1;
  std::vector<TH1D> againstElectronLooseMVA6_Tau1;
  std::vector<TH1D> againstElectronMediumMVA6_Tau1;
  std::vector<TH1D> againstElectronTightMVA6_Tau1;
  std::vector<TH1D> againstElectronVTightMVA6_Tau1;
  std::vector<TH1D> againstMuonLoose3_Tau1;
  std::vector<TH1D> againstMuonTight3_Tau1;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1;

  std::vector<TH1D> againstElectronVLooseMVA6_Tau2;
  std::vector<TH1D> againstElectronLooseMVA6_Tau2;
  std::vector<TH1D> againstElectronMediumMVA6_Tau2;
  std::vector<TH1D> againstElectronTightMVA6_Tau2;
  std::vector<TH1D> againstElectronVTightMVA6_Tau2;
  std::vector<TH1D> againstMuonLoose3_Tau2;
  std::vector<TH1D> againstMuonTight3_Tau2;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2;
  */
  std::vector<TH1D> ExtraLeptonVeto;
  std::vector<TH1D> Tau2HPSDecayMode;
  std::vector<TH1D> Tau1HPSDecayMode;

  std::vector<TH1D> TauTauVisMass;
  std::vector<TH1D> TauTauTruthMass;
  std::vector<TH1D> TauTauFullMass;
 
  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> TauTauFullMass_B;
  std::vector<TH1D> TauTauFullMass_C;
  std::vector<TH1D> TauTauFullMass_D;
  
  std::vector<TH1D> MET;
  std::vector<TH1D> METphi;
  std::vector<TH1D> PUPPImet;
  std::vector<TH1D> PUPPImetphi;
  std::vector<TH1D> TransverseMass;
  
  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> NPU;
  std::vector<TH1D> RHO;
  
  std::vector<TH1D> NbJets;
  
  std::vector<TH1D> h_SVFitMass; 
  std::vector<TH1D> h_SVFitStatus;
  /* std::vector<TH1D> svfTau1E; */
  /* std::vector<TH1D> svfTau2E; */
  
  /* std::vector<TH1D> PhiDatasvfitpipi; */
  /* std::vector<TH1D> PhiDatasvfitpirho; */
  /* std::vector<TH1D> PhiDatasvfitpia1; */
  /* std::vector<TH1D> PhiDatasvfitrhorho; */
  /* std::vector<TH1D> PhiDatasvfitrhoa1; */
  /* std::vector<TH1D> PhiDatasvfita1a1; */
  /*
  std::vector<TH1D> PhiDatavispipi;
  std::vector<TH1D> PhiDatavispirho;
  std::vector<TH1D> PhiDatavispia1;
  std::vector<TH1D> PhiDatavisrhorho;
  std::vector<TH1D> PhiDatavisrhoa1;
  std::vector<TH1D> PhiDatavisa1a1;
  */
  
  /* std::vector<TH1D> Etasvfit; */
  /* std::vector<TH1D> Phisvfitpipi; */
  /* std::vector<TH1D> Phisvfitpirho; */
  /* //std::vector<TH1D> Phisvfitlpi; */
  /* //std::vector<TH1D> Phisvfitlrho; */
  /* std::vector<TH1D> Phisvfitpia1; */
  /* std::vector<TH1D> Phisvfitrhorho; */
  /* std::vector<TH1D> Phisvfitrhoa1; */
  /* // std::vector<TH1D> Phisvfitla1; */
  /* std::vector<TH1D> Phisvfita1a1; */
  /* std::vector<TH1D> Thetasvfit; */
  
  std::vector<TH1D> Etavis;
  std::vector<TH1D> Phivispipi;
  std::vector<TH1D> Phivispirho;
  // std::vector<TH1D> Phivislpi;
  // std::vector<TH1D> Phivislrho;
  std::vector<TH1D> Phivispia1;
  std::vector<TH1D> Phivisrhorho;
  std::vector<TH1D> Phivisrhoa1;
  //std::vector<TH1D> Phivisla1;
  std::vector<TH1D> Phivisa1a1;
  std::vector<TH1D> Thetavis;
  
  std::vector<TH1D> Etatruth;
  std::vector<TH1D> Phitruthpipi;
  std::vector<TH1D> Phitruthpirho;
  // std::vector<TH1D> Phitruthlpi;
  // std::vector<TH1D> Phitruthlrho;
  std::vector<TH1D> Phitruthpia1;
  std::vector<TH1D> Phitruthrhorho;
  std::vector<TH1D> Phitruthrhoa1;
  // std::vector<TH1D> Phitruthla1;
  std::vector<TH1D> Phitrutha1a1;
  std::vector<TH1D> Thetatruth;
  /*
  std::vector<TH1D> PhiSvFitRespipi;
  std::vector<TH1D> PhiSvFitRespirho;
  std::vector<TH1D> PhiSvFitReslpi;
  std::vector<TH1D> PhiSvFitReslrho;
  std::vector<TH1D> PhiSvFitRespia1;
   std::vector<TH1D> PhiSvFitResrhorho;
  std::vector<TH1D> PhiSvFitResrhoa1;
  std::vector<TH1D> PhiSvFitResla1;
  std::vector<TH1D> PhiSvFitResa1a1;
  
  std::vector<TH1D> PhiVisRespipi;
  std::vector<TH1D> PhiVisRespirho;
  std::vector<TH1D> PhiVisReslpi;
  std::vector<TH1D> PhiVisReslrho;
  std::vector<TH1D> PhiVisRespia1;
  std::vector<TH1D> PhiVisResrhorho;
  std::vector<TH1D> PhiVisResrhoa1;
  std::vector<TH1D> PhiVisResla1;
  std::vector<TH1D> PhiVisResa1a1;  
  */

  /*
  std::vector<TH1D> TauTauFullPtRes;
  std::vector<TH1D> TauTauFullEtaRes;
  std::vector<TH1D> TauTauFullPhiRes;

  std::vector<TH1D> TauplusFullPtRes;
  std::vector<TH1D> TauplusFullEtaRes;
  std::vector<TH1D> TauplusFullPhiRes;

  std::vector<TH1D> TauminusFullPtRes;
  std::vector<TH1D> TauminusFullEtaRes;
  std::vector<TH1D> TauminusFullPhiRes;
  
  std::vector<TH1D> TauTauVisPtRes;
  std::vector<TH1D> TauTauVisEtaRes;
  std::vector<TH1D> TauTauVisPhiRes;

  std::vector<TH1D> TauplusVisPtRes;
  std::vector<TH1D> TauplusVisEtaRes;
  std::vector<TH1D> TauplusVisPhiRes;
    
  std::vector<TH1D> TauminusVisPtRes;
  std::vector<TH1D> TauminusVisEtaRes;
  std::vector<TH1D> TauminusVisPhiRes;

  std::vector<TH1D> TauTauFullPtResPull;
  std::vector<TH1D> TauTauFullEtaResPull;
  std::vector<TH1D> TauTauFullPhiResPull;

  std::vector<TH1D> TauplusFullPtResPull;
  std::vector<TH1D> TauplusFullEtaResPull;
  std::vector<TH1D> TauplusFullPhiResPull;

  std::vector<TH1D> TauminusFullPtResPull;
  std::vector<TH1D> TauminusFullEtaResPull;
  std::vector<TH1D> TauminusFullPhiResPull;

  std::vector<TH1D> TauTauVisPtResPull;
  std::vector<TH1D> TauTauVisEtaResPull;
  std::vector<TH1D> TauTauVisPhiResPull;

  std::vector<TH1D> TauplusVisPtResPull;
  std::vector<TH1D> TauplusVisEtaResPull;
  std::vector<TH1D> TauplusVisPhiResPull;
    
  std::vector<TH1D> TauminusVisPtResPull;
  std::vector<TH1D> TauminusVisEtaResPull;
  std::vector<TH1D> TauminusVisPhiResPull;
*/
  //std::vector<TH1D> DRTruth;
  //std::vector<TH1D> DRFull;
  // std::vector<TH1D> DRFullTruth;
  //std::vector<TH1D> DRVisTruth;
  
  std::vector<TH1D> Pi0EnergyRes;
  std::vector<TH1D> Pi0EnergyResPull;

  std::vector<TH1D> ZPtVis;

  std::vector<TH2D> NewPhivsDeltaPhi;
  std::vector<TH2D> NewPhivsDeltaEta;
  std::vector<TH2D> NewPhivsPhiproton;
  std::vector<TH2D> NewPhivsPhiTauplus;
  std::vector<TH2D> NewPhivsEtaproton;
  std::vector<TH2D> NewPhivsEtaTauplus;
  std::vector<TH2D> NewPhivsZPt;
  std::vector<TH1D> NewPhiSignal;
  std::vector<TH1D> NewPhiQCD;

};
#endif
