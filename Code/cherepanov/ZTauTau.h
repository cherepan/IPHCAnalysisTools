#ifndef ZTauTau_h
#define ZTauTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"
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
 private:
  // Selection Variables and Histos
  ClassicSVfit svfitAlgo1;
  // ClassicSVfit svfitAlgo2;
  //SVFitStorage svfitstorage;
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

  std::vector<TH1D> ExtraLeptonVeto;
  std::vector<TH1D> Tau2HPSDecayMode;
  std::vector<TH1D> Tau1HPSDecayMode;

  std::vector<TH1D> TauTauMass;

    std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;

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
  std::vector<TH1D> svfTau1E;
  std::vector<TH1D> svfTau2E;


};
#endif
