//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 15 17:49:09 2017 by ROOT version 5.34/18
// from TTree HTauTauTree/HTauTauTree
// found on file: /home-pbs/vcherepa/cms_work/CMSSW_8_0_25/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis.root
//////////////////////////////////////////////////////////

#ifndef NtupleReader_h
#define NtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   ULong64_t       EventNumber;
   Int_t           RunNumber;
   Int_t           lumi;
   Int_t           NBadMu;
   Long64_t        triggerbit;
   Int_t           metfilterbit;
   Float_t         met;
   Float_t         metphi;
   Float_t         PUPPImet;
   Float_t         PUPPImetphi;
   vector<float>   *daughters_IetaIeta;
   vector<float>   *daughters_full5x5_IetaIeta;
   vector<float>   *daughters_hOverE;
   vector<float>   *daughters_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *daughters_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *daughters_IoEmIoP;
   vector<float>   *daughters_IoEmIoP_ttH;
   Float_t         PFMETCov00;
   Float_t         PFMETCov01;
   Float_t         PFMETCov10;
   Float_t         PFMETCov11;
   Float_t         PFMETsignif;
   Int_t           npv;
   Float_t         npu;
   Float_t         PUReweight;
   Float_t         rho;
   vector<float>   *mothers_px;
   vector<float>   *mothers_py;
   vector<float>   *mothers_pz;
   vector<float>   *mothers_e;
   vector<Long64_t> *mothers_trgSeparateMatch;
   vector<string>  *trigger_name;
   vector<int>     *trigger_accept;
   vector<float>   *daughters_px;
   vector<float>   *daughters_py;
   vector<float>   *daughters_pz;
   vector<float>   *daughters_e;
   vector<int>     *daughters_charge;
   vector<float>   *daughters_charged_px;
   vector<float>   *daughters_charged_py;
   vector<float>   *daughters_charged_pz;
   vector<float>   *daughters_charged_e;
   vector<float>   *daughters_neutral_px;
   vector<float>   *daughters_neutral_py;
   vector<float>   *daughters_neutral_pz;
   vector<float>   *daughters_neutral_e;
   vector<int>     *daughters_TauUpExists;
   vector<float>   *daughters_px_TauUp;
   vector<float>   *daughters_py_TauUp;
   vector<float>   *daughters_pz_TauUp;
   vector<float>   *daughters_e_TauUp;
   vector<int>     *daughters_TauDownExists;
   vector<float>   *daughters_px_TauDown;
   vector<float>   *daughters_py_TauDown;
   vector<float>   *daughters_pz_TauDown;
   vector<float>   *daughters_e_TauDown;
   Int_t           PUNumInteractions;
   vector<int>     *daughters_genindex;
   Float_t         MC_weight;
   Float_t         MC_weight_scale_muF0p5;
   Float_t         MC_weight_scale_muF2;
   Float_t         MC_weight_scale_muR0p5;
   Float_t         MC_weight_scale_muR2;
   Float_t         lheHt;
   Int_t           lheNOutPartons;
   Int_t           lheNOutB;
   Int_t           lheNOutC;
   Float_t         aMCatNLOweight;
   vector<float>   *genpart_px;
   vector<float>   *genpart_py;
   vector<float>   *genpart_pz;
   vector<float>   *genpart_e;
   Int_t           DataMC_Type_idx;
   bool           Event_isRealData;
   vector<float>   *genpart_pca_x;
   vector<float>   *genpart_pca_y;
   vector<float>   *genpart_pca_z;
   vector<int>     *genpart_pdg;
   vector<int>     *genpart_status;
   vector<int>     *genpart_HMothInd;
   vector<int>     *genpart_MSSMHMothInd;
   vector<int>     *genpart_TopMothInd;
   vector<int>     *genpart_TauMothInd;
   vector<int>     *genpart_ZMothInd;
   vector<int>     *genpart_WMothInd;
   vector<int>     *genpart_bMothInd;
   vector<int>     *genpart_HZDecayMode;
   vector<int>     *genpart_TopDecayMode;
   vector<int>     *genpart_WDecayMode;
   vector<int>     *genpart_TauGenDecayMode;
   vector<int>     *genpart_TauGenDetailedDecayMode;
   vector<int>     *genpart_flags;
   vector<float>   *genjet_px;
   vector<float>   *genjet_py;
   vector<float>   *genjet_pz;
   vector<float>   *genjet_e;
   vector<int>     *genjet_partonFlavour;
   vector<int>     *genjet_hadronFlavour;
   Int_t           NUP;
   vector<float>   *SVfit_fitMETPhiTauUp;
   vector<float>   *SVfit_fitMETPhiTauDown;
   vector<float>   *SVfit_fitMETRhoTauUp;
   vector<float>   *SVfit_fitMETRhoTauDown;
   vector<float>   *SVfit_phiUncTauUp;
   vector<float>   *SVfit_phiUncTauDown;
   vector<float>   *SVfit_phiTauUp;
   vector<float>   *SVfit_phiTauDown;
   vector<float>   *SVfit_etaUncTauUp;
   vector<float>   *SVfit_etaUncTauDown;
   vector<float>   *SVfit_etaTauUp;
   vector<float>   *SVfit_etaTauDown;
   vector<float>   *SVfit_ptUncTauUp;
   vector<float>   *SVfit_ptUncTauDown;
   vector<float>   *SVfit_ptTauUp;
   vector<float>   *SVfit_ptTauDown;
   vector<float>   *SVfitTransverseMassTauUp;
   vector<float>   *SVfitTransverseMassTauDown;
   vector<float>   *SVfitMassTauUp;
   vector<float>   *SVfitMassTauDown;
   vector<vector<double> > *MCSignalParticle_p4;
   vector<int>     *MCSignalParticle_pdgid;
   vector<int>     *MCSignalParticle_charge;
   vector<vector<double> > *MCSignalParticle_Poca;
   vector<vector<unsigned int> > *MCSignalParticle_Tauidx;
   vector<vector<vector<double> > > *MCTauandProd_p4;
   vector<vector<vector<double> > > *MCTauandProd_Vertex;
   vector<vector<int> > *MCTauandProd_pdgid;
   vector<vector<unsigned int> > *MCTauandProd_midx;
   vector<vector<int> > *MCTauandProd_charge;
   vector<unsigned int> *MCTau_JAK;
   vector<unsigned int> *MCTau_DecayBitMask;
   vector<vector<float> > *MC_p4;
   vector<int>     *MC_pdgid;
   vector<int>     *MC_charge;
   vector<int>     *MC_midx;
   vector<vector<int> > *MC_childpdgid;
   vector<vector<int> > *MC_childidx;
   vector<int>     *MC_status;
   vector<float>   *SVfitMass;
   vector<float>   *SVfitTransverseMass;
   vector<float>   *SVfit_pt;
   vector<float>   *SVfit_ptUnc;
   vector<float>   *SVfit_eta;
   vector<float>   *SVfit_etaUnc;
   vector<float>   *SVfit_phi;
   vector<float>   *SVfit_phiUnc;
   vector<float>   *SVfit_fitMETRho;
   vector<float>   *SVfit_fitMETPhi;
   vector<bool>    *isOSCand;
   vector<float>   *METx;
   vector<float>   *METy;
   vector<float>   *uncorrMETx;
   vector<float>   *uncorrMETy;
   vector<float>   *MET_cov00;
   vector<float>   *MET_cov01;
   vector<float>   *MET_cov10;
   vector<float>   *MET_cov11;
   vector<float>   *MET_significance;
   vector<float>   *mT_Dau1;
   vector<float>   *mT_Dau2;
   vector<int>     *PDGIdDaughters;
   vector<int>     *indexDau1;
   vector<int>     *indexDau2;
   vector<int>     *particleType;
   vector<float>   *discriminator;
   vector<int>     *daughters_muonID;
   vector<int>     *daughters_typeOfMuon;
   vector<float>   *dxy;
   vector<float>   *dz;
   vector<float>   *dxy_innerTrack;
   vector<float>   *dz_innerTrack;
   vector<float>   *daughters_rel_error_trackpt;
   vector<float>   *SIP;
   vector<bool>    *daughters_iseleBDT;
   vector<bool>    *daughters_iseleWP80;
   vector<bool>    *daughters_iseleWP90;
   vector<float>   *daughters_eleMVAnt;
   vector<float>   *daughters_eleMVA_HZZ;
   vector<bool>    *daughters_passConversionVeto;
   vector<int>     *daughters_eleMissingHits;
   vector<bool>    *daughters_iseleChargeConsistent;
   vector<int>     *daughters_eleCUTID;
   vector<int>     *decayMode;
   vector<Long64_t> *tauID;
   vector<float>   *combreliso;
   vector<float>   *combreliso03;
   vector<float>   *daughters_depositR03_tracker;
   vector<float>   *daughters_depositR03_ecal;
   vector<float>   *daughters_depositR03_hcal;
   vector<int>     *daughters_decayModeFindingOldDMs;
   vector<float>   *daughters_SCeta;
   vector<float>   *againstElectronMVA5category;
   vector<float>   *againstElectronMVA5raw;
   vector<float>   *byPileupWeightedIsolationRaw3Hits;
   vector<float>   *footprintCorrection;
   vector<float>   *neutralIsoPtSumWeight;
   vector<float>   *photonPtSumOutsideSignalCone;
   vector<int>     *daughters_decayModeFindingNewDMs;
   vector<float>   *daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<float>   *daughters_byIsolationMVA3oldDMwoLTraw;
   vector<float>   *daughters_byIsolationMVA3oldDMwLTraw;
   vector<float>   *daughters_byIsolationMVA3newDMwoLTraw;
   vector<float>   *daughters_byIsolationMVA3newDMwLTraw;
   vector<float>   *daughters_byIsolationMVArun2v1DBoldDMwLTraw;
   vector<float>   *daughters_chargedIsoPtSum;
   vector<float>   *daughters_neutralIsoPtSum;
   vector<float>   *daughters_puCorrPtSum;
   vector<int>     *daughters_numChargedParticlesSignalCone;
   vector<int>     *daughters_numNeutralHadronsSignalCone;
   vector<int>     *daughters_numPhotonsSignalCone;
   vector<int>     *daughters_daughters_numParticlesSignalCone;
   vector<int>     *daughters_numChargedParticlesIsoCone;
   vector<int>     *daughters_numNeutralHadronsIsoCone;
   vector<int>     *daughters_numPhotonsIsoCone;
   vector<int>     *daughters_numParticlesIsoCone;
   vector<float>   *daughters_leadChargedParticlePt;
   vector<float>   *daughters_trackRefPt;
   vector<int>     *daughters_isLastTriggerObjectforPath;
   vector<Long64_t> *daughters_trgMatched;
   vector<int>     *daughters_isTriggerObjectforPath;
   vector<Long64_t> *daughters_FilterFired;
   vector<Long64_t> *daughters_isGoodTriggerType;
   vector<Long64_t> *daughters_L3FilterFired;
   vector<Long64_t> *daughters_L3FilterFiredLast;
   vector<float>   *daughters_HLTpt;
   vector<bool>    *daughters_isL1IsoTau28Matched;
   vector<int>     *Muon_trackCharge;
   vector<int>     *Muon_pdgid;
   vector<double>  *Muon_B;
   vector<double>  *Muon_M;
   vector<vector<double> > *Muon_par;
   vector<vector<double> > *Muon_cov;
   vector<vector<double> > *PFTauSVPos;
   vector<vector<double> > *PFTauSVCov;
   vector<vector<vector<double> > > *PFTauPionsP4;
   vector<vector<double> > *PFTauPionsCharge;
   vector<vector<double> > *PFTauSVChi2NDofMatchingQuality;
   vector<int>     *daughters_jetNDauChargedMVASel;
   vector<float>   *daughters_miniRelIsoCharged;
   vector<float>   *daughters_miniRelIsoNeutral;
   vector<float>   *daughters_jetPtRel;
   vector<float>   *daughters_jetPtRatio;
   vector<float>   *daughters_jetBTagCSV;
   vector<float>   *daughters_lepMVA_mvaId;
   vector<float>   *daughters_pca_x;
   vector<float>   *daughters_pca_y;
   vector<float>   *daughters_pca_z;
   vector<float>   *daughters_pcaRefitPV_x;
   vector<float>   *daughters_pcaRefitPV_y;
   vector<float>   *daughters_pcaRefitPV_z;
   vector<float>   *daughters_pcaGenPV_x;
   vector<float>   *daughters_pcaGenPV_y;
   vector<float>   *daughters_pcaGenPV_z;
   vector<vector<double> > *PFTau_a1_lvp;
   vector<vector<double> > *PFTau_a1_cov;
   vector<int>     *PFTau_a1_charge;
   vector<int>     *PFTau_a1_pdgid;
   vector<double>  *PFTau_a1_B;
   vector<double>  *PFTau_a1_M;
    Int_t           JetsNumber;
   vector<float>   *jets_px;
   vector<float>   *jets_py;
   vector<float>   *jets_pz;
   vector<float>   *jets_e;
   vector<float>   *jets_rawPt;
   vector<float>   *jets_area;
   vector<float>   *jets_mT;
   vector<int>     *jets_Flavour; 
   vector<int>     *jets_HadronFlavour;
   vector<int>     *jets_genjetIndex;
   vector<float>   *jets_PUJetID;
   vector<float>   *jets_PUJetIDupdated;
   vector<float>   *jets_vtxPt;
   vector<float>   *jets_vtxMass;
   vector<float>   *jets_vtx3dL;
   vector<float>   *jets_vtxNtrk;
   vector<float>   *jets_vtx3deL;
   vector<float>   *jets_leadTrackPt;
   vector<float>   *jets_leptonPtRel;
   vector<float>   *jets_leptonPt;
   vector<float>   *jets_leptonDeltaR;
   vector<float>   *jets_chEmEF;
   vector<float>   *jets_chHEF;
   vector<float>   *jets_nEmEF;
   vector<float>   *jets_nHEF;
   vector<float>   *jets_MUF;
   vector<int>     *jets_neMult;
   vector<int>     *jets_chMult;
   vector<float>   *jets_jecUnc;
   vector<float>   *bDiscriminator;
   vector<float>   *bCSVscore;
   vector<float>   *pfCombinedMVAV2BJetTags;
   vector<int>     *PFjetID;
   vector<float>   *jetRawf;
   vector<float>   *ak8jets_px;
   vector<float>   *ak8jets_py;
   vector<float>   *ak8jets_pz;
   vector<float>   *ak8jets_e;
   vector<float>   *ak8jets_SoftDropMass;
   vector<float>   *ak8jets_PrunedMass;
   vector<float>   *ak8jets_TrimmedMass;
   vector<float>   *ak8jets_FilteredMass;
   vector<float>   *ak8jets_tau1;
   vector<float>   *ak8jets_tau2;
   vector<float>   *ak8jets_tau3;
   vector<float>   *ak8jets_CSV;
   vector<int>     *ak8jets_nsubjets;
   vector<float>   *subjets_px;
   vector<float>   *subjets_py;
   vector<float>   *subjets_pz;
   vector<float>   *subjets_e;
   vector<float>   *subjets_CSV;
   vector<int>     *subjets_ak8MotherIdx;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   vector<double>  *pv_cov;
   Float_t         pvRefit_x;
   Float_t         pvRefit_y;
   Float_t         pvRefit_z;
   Float_t         pvGen_x;
   Float_t         pvGen_y;
   Float_t         pvGen_z;
   Bool_t          isRefitPV;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_NBadMu;   //!
   TBranch        *b_triggerbit;   //!
   TBranch        *b_metfilterbit;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_PUPPImet;   //!
   TBranch        *b_PUPPImetphi;   //!
   TBranch        *b_daughters_IetaIeta;   //!
   TBranch        *b_daughters_full5x5_IetaIeta;   //!
   TBranch        *b_daughters_hOverE;   //!
   TBranch        *b_daughters_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_daughters_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_daughters_IoEmIoP;   //!
   TBranch        *b_daughters_IoEmIoP_ttH;   //!
   TBranch        *b_PFMETCov00;   //!
   TBranch        *b_PFMETCov01;   //!
   TBranch        *b_PFMETCov10;   //!
   TBranch        *b_PFMETCov11;   //!
   TBranch        *b_PFMETsignif;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_PUReweight;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_mothers_px;   //!
   TBranch        *b_mothers_py;   //!
   TBranch        *b_mothers_pz;   //!
   TBranch        *b_mothers_e;   //!
   TBranch        *b_mothers_trgSeparateMatch;   //!
   TBranch        *b_trigger_name;   //!
   TBranch        *b_trigger_accept;   //!
   TBranch        *b_daughters_px;   //!
   TBranch        *b_daughters_py;   //!
   TBranch        *b_daughters_pz;   //!
   TBranch        *b_daughters_e;   //!
   TBranch        *b_daughters_charge;   //!
   TBranch        *b_daughters_charged_px;   //!
   TBranch        *b_daughters_charged_py;   //!
   TBranch        *b_daughters_charged_pz;   //!
   TBranch        *b_daughters_charged_e;   //!
   TBranch        *b_daughters_neutral_px;   //!
   TBranch        *b_daughters_neutral_py;   //!
   TBranch        *b_daughters_neutral_pz;   //!
   TBranch        *b_daughters_neutral_e;   //!
   TBranch        *b_daughters_TauUpExists;   //!
   TBranch        *b_daughters_px_TauUp;   //!
   TBranch        *b_daughters_py_TauUp;   //!
   TBranch        *b_daughters_pz_TauUp;   //!
   TBranch        *b_daughters_e_TauUp;   //!
   TBranch        *b_daughters_TauDownExists;   //!
   TBranch        *b_daughters_px_TauDown;   //!
   TBranch        *b_daughters_py_TauDown;   //!
   TBranch        *b_daughters_pz_TauDown;   //!
   TBranch        *b_daughters_e_TauDown;   //!
   TBranch        *b_PUNumInteractions;   //!
   TBranch        *b_daughters_genindex;   //!
   TBranch        *b_MC_weight;   //!
   TBranch        *b_MC_weight_scale_muF0p5;   //!
   TBranch        *b_MC_weight_scale_muF2;   //!
   TBranch        *b_MC_weight_scale_muR0p5;   //!
   TBranch        *b_MC_weight_scale_muR2;   //!
   TBranch        *b_lheHt;   //!
   TBranch        *b_lheNOutPartons;   //!
   TBranch        *b_lheNOutB;   //!
   TBranch        *b_lheNOutC;   //!
   TBranch        *b_aMCatNLOweight;   //!
   TBranch        *b_genpart_px;   //!
   TBranch        *b_genpart_py;   //!
   TBranch        *b_genpart_pz;   //!
   TBranch        *b_genpart_e;   //!
   TBranch        *b_DataMC_Type_idx;   //!
   TBranch        *b_Event_isRealData;   //!
   TBranch        *b_genpart_pca_x;   //!
   TBranch        *b_genpart_pca_y;   //!
   TBranch        *b_genpart_pca_z;   //!
   TBranch        *b_genpart_pdg;   //!
   TBranch        *b_genpart_status;   //!
   TBranch        *b_genpart_HMothInd;   //!
   TBranch        *b_genpart_MSSMHMothInd;   //!
   TBranch        *b_genpart_TopMothInd;   //!
   TBranch        *b_genpart_TauMothInd;   //!
   TBranch        *b_genpart_ZMothInd;   //!
   TBranch        *b_genpart_WMothInd;   //!
   TBranch        *b_genpart_bMothInd;   //!
   TBranch        *b_genpart_HZDecayMode;   //!
   TBranch        *b_genpart_TopDecayMode;   //!
   TBranch        *b_genpart_WDecayMode;   //!
   TBranch        *b_genpart_TauGenDecayMode;   //!
   TBranch        *b_genpart_TauGenDetailedDecayMode;   //!
   TBranch        *b_genpart_flags;   //!
   TBranch        *b_genjet_px;   //!
   TBranch        *b_genjet_py;   //!
   TBranch        *b_genjet_pz;   //!
   TBranch        *b_genjet_e;   //!
   TBranch        *b_genjet_partonFlavour;   //!
   TBranch        *b_genjet_hadronFlavour;   //!
   TBranch        *b_NUP;   //!
   TBranch        *b_SVfit_fitMETPhiTauUp;   //!
   TBranch        *b_SVfit_fitMETPhiTauDown;   //!
   TBranch        *b_SVfit_fitMETRhoTauUp;   //!
   TBranch        *b_SVfit_fitMETRhoTauDown;   //!
   TBranch        *b_SVfit_phiUncTauUp;   //!
   TBranch        *b_SVfit_phiUncTauDown;   //!
   TBranch        *b_SVfit_phiTauUp;   //!
   TBranch        *b_SVfit_phiTauDown;   //!
   TBranch        *b_SVfit_etaUncTauUp;   //!
   TBranch        *b_SVfit_etaUncTauDown;   //!
   TBranch        *b_SVfit_etaTauUp;   //!
   TBranch        *b_SVfit_etaTauDown;   //!
   TBranch        *b_SVfit_ptUncTauUp;   //!
   TBranch        *b_SVfit_ptUncTauDown;   //!
   TBranch        *b_SVfit_ptTauUp;   //!
   TBranch        *b_SVfit_ptTauDown;   //!
   TBranch        *b_SVfitTransverseMassTauUp;   //!
   TBranch        *b_SVfitTransverseMassTauDown;   //!
   TBranch        *b_SVfitMassTauUp;   //!
   TBranch        *b_SVfitMassTauDown;   //!
   TBranch        *b_MCSignalParticle_p4;   //!
   TBranch        *b_MCSignalParticle_pdgid;   //!
   TBranch        *b_MCSignalParticle_charge;   //!
   TBranch        *b_MCSignalParticle_Poca;   //!
   TBranch        *b_MCSignalParticle_Tauidx;   //!
   TBranch        *b_MCTauandProd_p4;   //!
   TBranch        *b_MCTauandProd_Vertex;   //!
   TBranch        *b_MCTauandProd_pdgid;   //!
   TBranch        *b_MCTauandProd_midx;   //!
   TBranch        *b_MCTauandProd_charge;   //!
   TBranch        *b_MCTau_JAK;   //!
   TBranch        *b_MCTau_DecayBitMask;   //!
   TBranch        *b_MC_p4;   //!
   TBranch        *b_MC_pdgid;   //!
   TBranch        *b_MC_charge;   //!
   TBranch        *b_MC_midx;   //!
   TBranch        *b_MC_childpdgid;   //!
   TBranch        *b_MC_childidx;   //!
   TBranch        *b_MC_status;   //!
   TBranch        *b_SVfitMass;   //!
   TBranch        *b_SVfitTransverseMass;   //!
   TBranch        *b_SVfit_pt;   //!
   TBranch        *b_SVfit_ptUnc;   //!
   TBranch        *b_SVfit_eta;   //!
   TBranch        *b_SVfit_etaUnc;   //!
   TBranch        *b_SVfit_phi;   //!
   TBranch        *b_SVfit_phiUnc;   //!
   TBranch        *b_SVfit_fitMETRho;   //!
   TBranch        *b_SVfit_fitMETPhi;   //!
   TBranch        *b_isOSCand;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_uncorrMETx;   //!
   TBranch        *b_uncorrMETy;   //!
   TBranch        *b_MET_cov00;   //!
   TBranch        *b_MET_cov01;   //!
   TBranch        *b_MET_cov10;   //!
   TBranch        *b_MET_cov11;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_mT_Dau1;   //!
   TBranch        *b_mT_Dau2;   //!
   TBranch        *b_PDGIdDaughters;   //!
   TBranch        *b_indexDau1;   //!
   TBranch        *b_indexDau2;   //!
   TBranch        *b_particleType;   //!
   TBranch        *b_discriminator;   //!
   TBranch        *b_daughters_muonID;   //!
   TBranch        *b_daughters_typeOfMuon;   //!
   TBranch        *b_dxy;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_dxy_innerTrack;   //!
   TBranch        *b_dz_innerTrack;   //!
   TBranch        *b_daughters_rel_error_trackpt;   //!
   TBranch        *b_SIP;   //!
   TBranch        *b_daughters_iseleBDT;   //!
   TBranch        *b_daughters_iseleWP80;   //!
   TBranch        *b_daughters_iseleWP90;   //!
   TBranch        *b_daughters_eleMVAnt;   //!
   TBranch        *b_daughters_eleMVA_HZZ;   //!
   TBranch        *b_daughters_passConversionVeto;   //!
   TBranch        *b_daughters_eleMissingHits;   //!
   TBranch        *b_daughters_iseleChargeConsistent;   //!
   TBranch        *b_daughters_eleCUTID;   //!
   TBranch        *b_decayMode;   //!
   TBranch        *b_tauID;   //!
   TBranch        *b_combreliso;   //!
   TBranch        *b_combreliso03;   //!
   TBranch        *b_daughters_depositR03_tracker;   //!
   TBranch        *b_daughters_depositR03_ecal;   //!
   TBranch        *b_daughters_depositR03_hcal;   //!
   TBranch        *b_daughters_decayModeFindingOldDMs;   //!
   TBranch        *b_daughters_SCeta;   //!
   TBranch        *b_againstElectronMVA5category;   //!
   TBranch        *b_againstElectronMVA5raw;   //!
   TBranch        *b_byPileupWeightedIsolationRaw3Hits;   //!
   TBranch        *b_footprintCorrection;   //!
   TBranch        *b_neutralIsoPtSumWeight;   //!
   TBranch        *b_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_daughters_decayModeFindingNewDMs;   //!
   TBranch        *b_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_daughters_byIsolationMVA3oldDMwoLTraw;   //!
   TBranch        *b_daughters_byIsolationMVA3oldDMwLTraw;   //!
   TBranch        *b_daughters_byIsolationMVA3newDMwoLTraw;   //!
   TBranch        *b_daughters_byIsolationMVA3newDMwLTraw;   //!
   TBranch        *b_daughters_byIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_daughters_chargedIsoPtSum;   //!
   TBranch        *b_daughters_neutralIsoPtSum;   //!
   TBranch        *b_daughters_puCorrPtSum;   //!
   TBranch        *b_daughters_numChargedParticlesSignalCone;   //!
   TBranch        *b_daughters_numNeutralHadronsSignalCone;   //!
   TBranch        *b_daughters_numPhotonsSignalCone;   //!
   TBranch        *b_daughters_daughters_numParticlesSignalCone;   //!
   TBranch        *b_daughters_numChargedParticlesIsoCone;   //!
   TBranch        *b_daughters_numNeutralHadronsIsoCone;   //!
   TBranch        *b_daughters_numPhotonsIsoCone;   //!
   TBranch        *b_daughters_numParticlesIsoCone;   //!
   TBranch        *b_daughters_leadChargedParticlePt;   //!
   TBranch        *b_daughters_trackRefPt;   //!
   TBranch        *b_daughters_isLastTriggerObjectforPath;   //!
   TBranch        *b_daughters_trgMatched;   //!
   TBranch        *b_daughters_isTriggerObjectforPath;   //!
   TBranch        *b_daughters_FilterFired;   //!
   TBranch        *b_daughters_isGoodTriggerType;   //!
   TBranch        *b_daughters_L3FilterFired;   //!
   TBranch        *b_daughters_L3FilterFiredLast;   //!
   TBranch        *b_daughters_HLTpt;   //!
   TBranch        *b_daughters_isL1IsoTau28Matched;   //!
   TBranch        *b_Muon_trackCharge;   //!
   TBranch        *b_Muon_pdgid;   //!
   TBranch        *b_Muon_B;   //!
   TBranch        *b_Muon_M;   //!
   TBranch        *b_Muon_par;   //!
   TBranch        *b_Muon_cov;   //!
    TBranch        *b_PFTauSVPos;   //!
   TBranch        *b_PFTauSVCov;   //!
   TBranch        *b_PFTauPionsP4;   //!
   TBranch        *b_PFTauPionsCharge;   //!
   TBranch        *b_PFTauSVChi2NDofMatchingQuality;   //!
   TBranch        *b_daughters_jetNDauChargedMVASel;   //!
   TBranch        *b_daughters_miniRelIsoCharged;   //!
   TBranch        *b_daughters_miniRelIsoNeutral;   //!
   TBranch        *b_daughters_jetPtRel;   //!
   TBranch        *b_daughters_jetPtRatio;   //!
   TBranch        *b_daughters_jetBTagCSV;   //!
   TBranch        *b_daughters_lepMVA_mvaId;   //!
   TBranch        *b_daughters_pca_x;   //!
   TBranch        *b_daughters_pca_y;   //!
   TBranch        *b_daughters_pca_z;   //!
   TBranch        *b_daughters_pcaRefitPV_x;   //!
   TBranch        *b_daughters_pcaRefitPV_y;   //!
   TBranch        *b_daughters_pcaRefitPV_z;   //!
   TBranch        *b_daughters_pcaGenPV_x;   //!
   TBranch        *b_daughters_pcaGenPV_y;   //!
   TBranch        *b_daughters_pcaGenPV_z;   //!
   TBranch        *b_PFTau_a1_lvp;   //!
   TBranch        *b_PFTau_a1_cov;   //!
   TBranch        *b_PFTau_a1_charge;   //!
   TBranch        *b_PFTau_a1_pdgid;   //!
   TBranch        *b_PFTau_a1_B;   //!
   TBranch        *b_PFTau_a1_M;   //!
   TBranch        *b_JetsNumber;   //!
   TBranch        *b_jets_px;   //!
   TBranch        *b_jets_py;   //!
   TBranch        *b_jets_pz;   //!
   TBranch        *b_jets_e;   //!
   TBranch        *b_jets_rawPt;   //!
   TBranch        *b_jets_area;   //!
   TBranch        *b_jets_mT;   //!
   TBranch        *b_jets_Flavour;   //!
   TBranch        *b_jets_HadronFlavour;   //!
   TBranch        *b_jets_genjetIndex;   //!
   TBranch        *b_jets_PUJetID;   //!
   TBranch        *b_jets_PUJetIDupdated;   //!
   TBranch        *b_jets_vtxPt;   //!
   TBranch        *b_jets_vtxMass;   //!
   TBranch        *b_jets_vtx3dL;   //!
   TBranch        *b_jets_vtxNtrk;   //!
   TBranch        *b_jets_vtx3deL;   //!
   TBranch        *b_jets_leadTrackPt;   //!
   TBranch        *b_jets_leptonPtRel;   //!
   TBranch        *b_jets_leptonPt;   //!
   TBranch        *b_jets_leptonDeltaR;   //!
   TBranch        *b_jets_chEmEF;   //!
   TBranch        *b_jets_chHEF;   //!
   TBranch        *b_jets_nEmEF;   //!
   TBranch        *b_jets_nHEF;   //!
   TBranch        *b_jets_MUF;   //!
   TBranch        *b_jets_neMult;   //!
   TBranch        *b_jets_chMult;   //!
   TBranch        *b_jets_jecUnc;   //!
   TBranch        *b_bDiscriminator;   //!
   TBranch        *b_bCSVscore;   //!
   TBranch        *b_pfCombinedMVAV2BJetTags;   //!
   TBranch        *b_PFjetID;   //!
   TBranch        *b_jetRawf;   //!
   TBranch        *b_ak8jets_px;   //!
   TBranch        *b_ak8jets_py;   //!
   TBranch        *b_ak8jets_pz;   //!
   TBranch        *b_ak8jets_e;   //!
   TBranch        *b_ak8jets_SoftDropMass;   //!
   TBranch        *b_ak8jets_PrunedMass;   //!
   TBranch        *b_ak8jets_TrimmedMass;   //!
   TBranch        *b_ak8jets_FilteredMass;   //!
   TBranch        *b_ak8jets_tau1;   //!
   TBranch        *b_ak8jets_tau2;   //!
   TBranch        *b_ak8jets_tau3;   //!
   TBranch        *b_ak8jets_CSV;   //!
   TBranch        *b_ak8jets_nsubjets;   //!
   TBranch        *b_subjets_px;   //!
   TBranch        *b_subjets_py;   //!
   TBranch        *b_subjets_pz;   //!
   TBranch        *b_subjets_e;   //!
   TBranch        *b_subjets_CSV;   //!
   TBranch        *b_subjets_ak8MotherIdx;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_cov;   //!
   TBranch        *b_pvRefit_x;   //!
   TBranch        *b_pvRefit_y;   //!
   TBranch        *b_pvRefit_z;   //!
   TBranch        *b_pvGen_x;   //!
   TBranch        *b_pvGen_y;   //!
   TBranch        *b_pvGen_z;   //!
   TBranch        *b_isRefitPV;   //!

   NtupleReader(TTree *tree=0);
   virtual ~NtupleReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleReader_cxx
NtupleReader::NtupleReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {


#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("HTauTauTree/HTauTauTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("HTauTauTree/HTauTauTree","");
      chain->Add("/home-pbs/vcherepa/cms_work/CMSSW_8_0_25/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis.root/HTauTauTree/HTauTauTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}



NtupleReader::~NtupleReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupleReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupleReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtupleReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   daughters_IetaIeta = 0;
   daughters_full5x5_IetaIeta = 0;
   daughters_hOverE = 0;
   daughters_deltaEtaSuperClusterTrackAtVtx = 0;
   daughters_deltaPhiSuperClusterTrackAtVtx = 0;
   daughters_IoEmIoP = 0;
   daughters_IoEmIoP_ttH = 0;
   mothers_px = 0;
   mothers_py = 0;
   mothers_pz = 0;
   mothers_e = 0;
   mothers_trgSeparateMatch = 0;
   trigger_name = 0;
   trigger_accept = 0;
   daughters_px = 0;
   daughters_py = 0;
   daughters_pz = 0;
   daughters_e = 0;
   daughters_charge = 0;
   daughters_charged_px = 0;
   daughters_charged_py = 0;
   daughters_charged_pz = 0;
   daughters_charged_e = 0;
   daughters_neutral_px = 0;
   daughters_neutral_py = 0;
   daughters_neutral_pz = 0;
   daughters_neutral_e = 0;
   daughters_TauUpExists = 0;
   daughters_px_TauUp = 0;
   daughters_py_TauUp = 0;
   daughters_pz_TauUp = 0;
   daughters_e_TauUp = 0;
   daughters_TauDownExists = 0;
   daughters_px_TauDown = 0;
   daughters_py_TauDown = 0;
   daughters_pz_TauDown = 0;
   daughters_e_TauDown = 0;
   daughters_genindex = 0;
   genpart_px = 0;
   genpart_py = 0;
   genpart_pz = 0;
   genpart_e = 0;
   genpart_pca_x = 0;
   genpart_pca_y = 0;
   genpart_pca_z = 0;
   genpart_pdg = 0;
   genpart_status = 0;
   genpart_HMothInd = 0;
   genpart_MSSMHMothInd = 0;
   genpart_TopMothInd = 0;
   genpart_TauMothInd = 0;
   genpart_ZMothInd = 0;
   genpart_WMothInd = 0;
   genpart_bMothInd = 0;
   genpart_HZDecayMode = 0;
   genpart_TopDecayMode = 0;
   genpart_WDecayMode = 0;
   genpart_TauGenDecayMode = 0;
   genpart_TauGenDetailedDecayMode = 0;
   genpart_flags = 0;
   genjet_px = 0;
   genjet_py = 0;
   genjet_pz = 0;
   genjet_e = 0;
   genjet_partonFlavour = 0;
   genjet_hadronFlavour = 0;
   SVfit_fitMETPhiTauUp = 0;
   SVfit_fitMETPhiTauDown = 0;
   SVfit_fitMETRhoTauUp = 0;
   SVfit_fitMETRhoTauDown = 0;
   SVfit_phiUncTauUp = 0;
   SVfit_phiUncTauDown = 0;
   SVfit_phiTauUp = 0;
   SVfit_phiTauDown = 0;
   SVfit_etaUncTauUp = 0;
   SVfit_etaUncTauDown = 0;
   SVfit_etaTauUp = 0;
   SVfit_etaTauDown = 0;
   SVfit_ptUncTauUp = 0;
   SVfit_ptUncTauDown = 0;
   SVfit_ptTauUp = 0;
   SVfit_ptTauDown = 0;
   SVfitTransverseMassTauUp = 0;
   SVfitTransverseMassTauDown = 0;
   SVfitMassTauUp = 0;
   SVfitMassTauDown = 0;
   MCSignalParticle_p4 = 0;
   MCSignalParticle_pdgid = 0;
   MCSignalParticle_charge = 0;
   MCSignalParticle_Poca = 0;
   MCSignalParticle_Tauidx = 0;
   MCTauandProd_p4 = 0;
   MCTauandProd_Vertex = 0;
   MCTauandProd_pdgid = 0;
   MCTauandProd_midx = 0;
   MCTauandProd_charge = 0;
   MCTau_JAK = 0;
   MCTau_DecayBitMask = 0;
   MC_p4 = 0;
   MC_pdgid = 0;
   MC_charge = 0;
   MC_midx = 0;
   MC_childpdgid = 0;
   MC_childidx = 0;
   MC_status = 0;
   SVfitMass = 0;
   SVfitTransverseMass = 0;
   SVfit_pt = 0;
   SVfit_ptUnc = 0;
   SVfit_eta = 0;
   SVfit_etaUnc = 0;
   SVfit_phi = 0;
   SVfit_phiUnc = 0;
   SVfit_fitMETRho = 0;
   SVfit_fitMETPhi = 0;
   isOSCand = 0;
   METx = 0;
   METy = 0;
   uncorrMETx = 0;
   uncorrMETy = 0;
   MET_cov00 = 0;
   MET_cov01 = 0;
   MET_cov10 = 0;
   MET_cov11 = 0;
   MET_significance = 0;
   mT_Dau1 = 0;
   mT_Dau2 = 0;
   PDGIdDaughters = 0;
   indexDau1 = 0;
   indexDau2 = 0;
   particleType = 0;
   discriminator = 0;
   daughters_muonID = 0;
   daughters_typeOfMuon = 0;
   dxy = 0;
   dz = 0;
   dxy_innerTrack = 0;
   dz_innerTrack = 0;
   daughters_rel_error_trackpt = 0;
   SIP = 0;
   daughters_iseleBDT = 0;
   daughters_iseleWP80 = 0;
   daughters_iseleWP90 = 0;
   daughters_eleMVAnt = 0;
   daughters_eleMVA_HZZ = 0;
   daughters_passConversionVeto = 0;
   daughters_eleMissingHits = 0;
   daughters_iseleChargeConsistent = 0;
   daughters_eleCUTID = 0;
   decayMode = 0;
   tauID = 0;
   combreliso = 0;
   combreliso03 = 0;
   daughters_depositR03_tracker = 0;
   daughters_depositR03_ecal = 0;
   daughters_depositR03_hcal = 0;
   daughters_decayModeFindingOldDMs = 0;
   daughters_SCeta = 0;
   againstElectronMVA5category = 0;
   againstElectronMVA5raw = 0;
   byPileupWeightedIsolationRaw3Hits = 0;
   footprintCorrection = 0;
   neutralIsoPtSumWeight = 0;
   photonPtSumOutsideSignalCone = 0;
   daughters_decayModeFindingNewDMs = 0;
   daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   daughters_byIsolationMVA3oldDMwoLTraw = 0;
   daughters_byIsolationMVA3oldDMwLTraw = 0;
   daughters_byIsolationMVA3newDMwoLTraw = 0;
   daughters_byIsolationMVA3newDMwLTraw = 0;
   daughters_byIsolationMVArun2v1DBoldDMwLTraw = 0;
   daughters_chargedIsoPtSum = 0;
   daughters_neutralIsoPtSum = 0;
   daughters_puCorrPtSum = 0;
   daughters_numChargedParticlesSignalCone = 0;
   daughters_numNeutralHadronsSignalCone = 0;
   daughters_numPhotonsSignalCone = 0;
   daughters_daughters_numParticlesSignalCone = 0;
   daughters_numChargedParticlesIsoCone = 0;
   daughters_numNeutralHadronsIsoCone = 0;
   daughters_numPhotonsIsoCone = 0;
   daughters_numParticlesIsoCone = 0;
   daughters_leadChargedParticlePt = 0;
   daughters_trackRefPt = 0;
   daughters_isLastTriggerObjectforPath = 0;
   daughters_trgMatched = 0;
   daughters_isTriggerObjectforPath = 0;
   daughters_FilterFired = 0;
   daughters_isGoodTriggerType = 0;
   daughters_L3FilterFired = 0;
   daughters_L3FilterFiredLast = 0;
   daughters_HLTpt = 0;
   daughters_isL1IsoTau28Matched = 0;
   Muon_trackCharge = 0;
   Muon_pdgid = 0;
   Muon_B = 0;
   Muon_M = 0;
   Muon_par = 0;
   Muon_cov = 0;
   PFTauSVPos = 0;
   PFTauSVCov = 0;
   PFTauPionsP4 = 0;
   PFTauPionsCharge = 0;
   PFTauSVChi2NDofMatchingQuality = 0;
   daughters_jetNDauChargedMVASel = 0;
   daughters_miniRelIsoCharged = 0;
   daughters_miniRelIsoNeutral = 0;
   daughters_jetPtRel = 0;
   daughters_jetPtRatio = 0;
   daughters_jetBTagCSV = 0;
   daughters_lepMVA_mvaId = 0;
   daughters_pca_x = 0;
   daughters_pca_y = 0;
   daughters_pca_z = 0;
   daughters_pcaRefitPV_x = 0;
   daughters_pcaRefitPV_y = 0;
   daughters_pcaRefitPV_z = 0;
   daughters_pcaGenPV_x = 0;
   daughters_pcaGenPV_y = 0;
   daughters_pcaGenPV_z = 0;
   PFTau_a1_lvp = 0;
   PFTau_a1_cov = 0;
   PFTau_a1_charge = 0;
   PFTau_a1_pdgid = 0;
   PFTau_a1_B = 0;
   PFTau_a1_M = 0;
   JetsNumber = 0;
   jets_px = 0;
   jets_py = 0;
   jets_pz = 0;
   jets_e = 0;
   jets_rawPt = 0;
   jets_area = 0;
   jets_mT = 0;
   jets_Flavour = 0;
   jets_HadronFlavour = 0;
   jets_genjetIndex = 0;
   jets_PUJetID = 0;
   jets_PUJetIDupdated = 0;
   jets_vtxPt = 0;
   jets_vtxMass = 0;
   jets_vtx3dL = 0;
   jets_vtxNtrk = 0;
   jets_vtx3deL = 0;
   jets_leadTrackPt = 0;
   jets_leptonPtRel = 0;
   jets_leptonPt = 0;
   jets_leptonDeltaR = 0;
   jets_chEmEF = 0;
   jets_chHEF = 0;
   jets_nEmEF = 0;
   jets_nHEF = 0;
   jets_MUF = 0;
   jets_neMult = 0;
   jets_chMult = 0;
   jets_jecUnc = 0;
   bDiscriminator = 0;
   bCSVscore = 0;
   pfCombinedMVAV2BJetTags = 0;
   PFjetID = 0;
   jetRawf = 0;
   ak8jets_px = 0;
   ak8jets_py = 0;
   ak8jets_pz = 0;
   ak8jets_e = 0;
   ak8jets_SoftDropMass = 0;
   ak8jets_PrunedMass = 0;
   ak8jets_TrimmedMass = 0;
   ak8jets_FilteredMass = 0;
   ak8jets_tau1 = 0;
   ak8jets_tau2 = 0;
   ak8jets_tau3 = 0;
   ak8jets_CSV = 0;
   ak8jets_nsubjets = 0;
   subjets_px = 0;
   subjets_py = 0;
   subjets_pz = 0;
   subjets_e = 0;
   subjets_CSV = 0;
   subjets_ak8MotherIdx = 0;
   pv_cov = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("NBadMu", &NBadMu, &b_NBadMu);
   fChain->SetBranchAddress("triggerbit", &triggerbit, &b_triggerbit);
   fChain->SetBranchAddress("metfilterbit", &metfilterbit, &b_metfilterbit);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("PUPPImet", &PUPPImet, &b_PUPPImet);
   fChain->SetBranchAddress("PUPPImetphi", &PUPPImetphi, &b_PUPPImetphi);
   fChain->SetBranchAddress("daughters_IetaIeta", &daughters_IetaIeta, &b_daughters_IetaIeta);
   fChain->SetBranchAddress("daughters_full5x5_IetaIeta", &daughters_full5x5_IetaIeta, &b_daughters_full5x5_IetaIeta);
   fChain->SetBranchAddress("daughters_hOverE", &daughters_hOverE, &b_daughters_hOverE);
   fChain->SetBranchAddress("daughters_deltaEtaSuperClusterTrackAtVtx", &daughters_deltaEtaSuperClusterTrackAtVtx, &b_daughters_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("daughters_deltaPhiSuperClusterTrackAtVtx", &daughters_deltaPhiSuperClusterTrackAtVtx, &b_daughters_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("daughters_IoEmIoP", &daughters_IoEmIoP, &b_daughters_IoEmIoP);
   fChain->SetBranchAddress("daughters_IoEmIoP_ttH", &daughters_IoEmIoP_ttH, &b_daughters_IoEmIoP_ttH);
   fChain->SetBranchAddress("PFMETCov00", &PFMETCov00, &b_PFMETCov00);
   fChain->SetBranchAddress("PFMETCov01", &PFMETCov01, &b_PFMETCov01);
   fChain->SetBranchAddress("PFMETCov10", &PFMETCov10, &b_PFMETCov10);
   fChain->SetBranchAddress("PFMETCov11", &PFMETCov11, &b_PFMETCov11);
   fChain->SetBranchAddress("PFMETsignif", &PFMETsignif, &b_PFMETsignif);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("PUReweight", &PUReweight, &b_PUReweight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("mothers_px", &mothers_px, &b_mothers_px);
   fChain->SetBranchAddress("mothers_py", &mothers_py, &b_mothers_py);
   fChain->SetBranchAddress("mothers_pz", &mothers_pz, &b_mothers_pz);
   fChain->SetBranchAddress("mothers_e", &mothers_e, &b_mothers_e);
   fChain->SetBranchAddress("mothers_trgSeparateMatch", &mothers_trgSeparateMatch, &b_mothers_trgSeparateMatch);
   fChain->SetBranchAddress("trigger_name", &trigger_name, &b_trigger_name);
   fChain->SetBranchAddress("trigger_accept", &trigger_accept, &b_trigger_accept);
   fChain->SetBranchAddress("daughters_px", &daughters_px, &b_daughters_px);
   fChain->SetBranchAddress("daughters_py", &daughters_py, &b_daughters_py);
   fChain->SetBranchAddress("daughters_pz", &daughters_pz, &b_daughters_pz);
   fChain->SetBranchAddress("daughters_e", &daughters_e, &b_daughters_e);
   fChain->SetBranchAddress("daughters_charge", &daughters_charge, &b_daughters_charge);
   fChain->SetBranchAddress("daughters_charged_px", &daughters_charged_px, &b_daughters_charged_px);
   fChain->SetBranchAddress("daughters_charged_py", &daughters_charged_py, &b_daughters_charged_py);
   fChain->SetBranchAddress("daughters_charged_pz", &daughters_charged_pz, &b_daughters_charged_pz);
   fChain->SetBranchAddress("daughters_charged_e", &daughters_charged_e, &b_daughters_charged_e);
   fChain->SetBranchAddress("daughters_neutral_px", &daughters_neutral_px, &b_daughters_neutral_px);
   fChain->SetBranchAddress("daughters_neutral_py", &daughters_neutral_py, &b_daughters_neutral_py);
   fChain->SetBranchAddress("daughters_neutral_pz", &daughters_neutral_pz, &b_daughters_neutral_pz);
   fChain->SetBranchAddress("daughters_neutral_e", &daughters_neutral_e, &b_daughters_neutral_e);
   fChain->SetBranchAddress("daughters_TauUpExists", &daughters_TauUpExists, &b_daughters_TauUpExists);
   fChain->SetBranchAddress("daughters_px_TauUp", &daughters_px_TauUp, &b_daughters_px_TauUp);
   fChain->SetBranchAddress("daughters_py_TauUp", &daughters_py_TauUp, &b_daughters_py_TauUp);
   fChain->SetBranchAddress("daughters_pz_TauUp", &daughters_pz_TauUp, &b_daughters_pz_TauUp);
   fChain->SetBranchAddress("daughters_e_TauUp", &daughters_e_TauUp, &b_daughters_e_TauUp);
   fChain->SetBranchAddress("daughters_TauDownExists", &daughters_TauDownExists, &b_daughters_TauDownExists);
   fChain->SetBranchAddress("daughters_px_TauDown", &daughters_px_TauDown, &b_daughters_px_TauDown);
   fChain->SetBranchAddress("daughters_py_TauDown", &daughters_py_TauDown, &b_daughters_py_TauDown);
   fChain->SetBranchAddress("daughters_pz_TauDown", &daughters_pz_TauDown, &b_daughters_pz_TauDown);
   fChain->SetBranchAddress("daughters_e_TauDown", &daughters_e_TauDown, &b_daughters_e_TauDown);
   fChain->SetBranchAddress("PUNumInteractions", &PUNumInteractions, &b_PUNumInteractions);
   fChain->SetBranchAddress("daughters_genindex", &daughters_genindex, &b_daughters_genindex);
   fChain->SetBranchAddress("MC_weight", &MC_weight, &b_MC_weight);
   fChain->SetBranchAddress("MC_weight_scale_muF0p5", &MC_weight_scale_muF0p5, &b_MC_weight_scale_muF0p5);
   fChain->SetBranchAddress("MC_weight_scale_muF2", &MC_weight_scale_muF2, &b_MC_weight_scale_muF2);
   fChain->SetBranchAddress("MC_weight_scale_muR0p5", &MC_weight_scale_muR0p5, &b_MC_weight_scale_muR0p5);
   fChain->SetBranchAddress("MC_weight_scale_muR2", &MC_weight_scale_muR2, &b_MC_weight_scale_muR2);
   fChain->SetBranchAddress("lheHt", &lheHt, &b_lheHt);
   fChain->SetBranchAddress("lheNOutPartons", &lheNOutPartons, &b_lheNOutPartons);
   fChain->SetBranchAddress("lheNOutB", &lheNOutB, &b_lheNOutB);
   fChain->SetBranchAddress("lheNOutC", &lheNOutC, &b_lheNOutC);
   fChain->SetBranchAddress("aMCatNLOweight", &aMCatNLOweight, &b_aMCatNLOweight);
   fChain->SetBranchAddress("genpart_px", &genpart_px, &b_genpart_px);
   fChain->SetBranchAddress("genpart_py", &genpart_py, &b_genpart_py);
   fChain->SetBranchAddress("genpart_pz", &genpart_pz, &b_genpart_pz);
   fChain->SetBranchAddress("genpart_e", &genpart_e, &b_genpart_e);
   fChain->SetBranchAddress("DataMC_Type_idx", &DataMC_Type_idx, &b_DataMC_Type_idx);
   fChain->SetBranchAddress("Event_isRealData", &Event_isRealData, &b_Event_isRealData);
   fChain->SetBranchAddress("genpart_pca_x", &genpart_pca_x, &b_genpart_pca_x);
   fChain->SetBranchAddress("genpart_pca_y", &genpart_pca_y, &b_genpart_pca_y);
   fChain->SetBranchAddress("genpart_pca_z", &genpart_pca_z, &b_genpart_pca_z);
   fChain->SetBranchAddress("genpart_pdg", &genpart_pdg, &b_genpart_pdg);
   fChain->SetBranchAddress("genpart_status", &genpart_status, &b_genpart_status);
   fChain->SetBranchAddress("genpart_HMothInd", &genpart_HMothInd, &b_genpart_HMothInd);
   fChain->SetBranchAddress("genpart_MSSMHMothInd", &genpart_MSSMHMothInd, &b_genpart_MSSMHMothInd);
   fChain->SetBranchAddress("genpart_TopMothInd", &genpart_TopMothInd, &b_genpart_TopMothInd);
   fChain->SetBranchAddress("genpart_TauMothInd", &genpart_TauMothInd, &b_genpart_TauMothInd);
   fChain->SetBranchAddress("genpart_ZMothInd", &genpart_ZMothInd, &b_genpart_ZMothInd);
   fChain->SetBranchAddress("genpart_WMothInd", &genpart_WMothInd, &b_genpart_WMothInd);
   fChain->SetBranchAddress("genpart_bMothInd", &genpart_bMothInd, &b_genpart_bMothInd);
   fChain->SetBranchAddress("genpart_HZDecayMode", &genpart_HZDecayMode, &b_genpart_HZDecayMode);
   fChain->SetBranchAddress("genpart_TopDecayMode", &genpart_TopDecayMode, &b_genpart_TopDecayMode);
   fChain->SetBranchAddress("genpart_WDecayMode", &genpart_WDecayMode, &b_genpart_WDecayMode);
   fChain->SetBranchAddress("genpart_TauGenDecayMode", &genpart_TauGenDecayMode, &b_genpart_TauGenDecayMode);
   fChain->SetBranchAddress("genpart_TauGenDetailedDecayMode", &genpart_TauGenDetailedDecayMode, &b_genpart_TauGenDetailedDecayMode);
   fChain->SetBranchAddress("genpart_flags", &genpart_flags, &b_genpart_flags);
   fChain->SetBranchAddress("genjet_px", &genjet_px, &b_genjet_px);
   fChain->SetBranchAddress("genjet_py", &genjet_py, &b_genjet_py);
   fChain->SetBranchAddress("genjet_pz", &genjet_pz, &b_genjet_pz);
   fChain->SetBranchAddress("genjet_e", &genjet_e, &b_genjet_e);
   fChain->SetBranchAddress("genjet_partonFlavour", &genjet_partonFlavour, &b_genjet_partonFlavour);
   fChain->SetBranchAddress("genjet_hadronFlavour", &genjet_hadronFlavour, &b_genjet_hadronFlavour);
   fChain->SetBranchAddress("NUP", &NUP, &b_NUP);
   fChain->SetBranchAddress("SVfit_fitMETPhiTauUp", &SVfit_fitMETPhiTauUp, &b_SVfit_fitMETPhiTauUp);
   fChain->SetBranchAddress("SVfit_fitMETPhiTauDown", &SVfit_fitMETPhiTauDown, &b_SVfit_fitMETPhiTauDown);
   fChain->SetBranchAddress("SVfit_fitMETRhoTauUp", &SVfit_fitMETRhoTauUp, &b_SVfit_fitMETRhoTauUp);
   fChain->SetBranchAddress("SVfit_fitMETRhoTauDown", &SVfit_fitMETRhoTauDown, &b_SVfit_fitMETRhoTauDown);
   fChain->SetBranchAddress("SVfit_phiUncTauUp", &SVfit_phiUncTauUp, &b_SVfit_phiUncTauUp);
   fChain->SetBranchAddress("SVfit_phiUncTauDown", &SVfit_phiUncTauDown, &b_SVfit_phiUncTauDown);
   fChain->SetBranchAddress("SVfit_phiTauUp", &SVfit_phiTauUp, &b_SVfit_phiTauUp);
   fChain->SetBranchAddress("SVfit_phiTauDown", &SVfit_phiTauDown, &b_SVfit_phiTauDown);
   fChain->SetBranchAddress("SVfit_etaUncTauUp", &SVfit_etaUncTauUp, &b_SVfit_etaUncTauUp);
   fChain->SetBranchAddress("SVfit_etaUncTauDown", &SVfit_etaUncTauDown, &b_SVfit_etaUncTauDown);
   fChain->SetBranchAddress("SVfit_etaTauUp", &SVfit_etaTauUp, &b_SVfit_etaTauUp);
   fChain->SetBranchAddress("SVfit_etaTauDown", &SVfit_etaTauDown, &b_SVfit_etaTauDown);
   fChain->SetBranchAddress("SVfit_ptUncTauUp", &SVfit_ptUncTauUp, &b_SVfit_ptUncTauUp);
   fChain->SetBranchAddress("SVfit_ptUncTauDown", &SVfit_ptUncTauDown, &b_SVfit_ptUncTauDown);
   fChain->SetBranchAddress("SVfit_ptTauUp", &SVfit_ptTauUp, &b_SVfit_ptTauUp);
   fChain->SetBranchAddress("SVfit_ptTauDown", &SVfit_ptTauDown, &b_SVfit_ptTauDown);
   fChain->SetBranchAddress("SVfitTransverseMassTauUp", &SVfitTransverseMassTauUp, &b_SVfitTransverseMassTauUp);
   fChain->SetBranchAddress("SVfitTransverseMassTauDown", &SVfitTransverseMassTauDown, &b_SVfitTransverseMassTauDown);
   fChain->SetBranchAddress("SVfitMassTauUp", &SVfitMassTauUp, &b_SVfitMassTauUp);
   fChain->SetBranchAddress("SVfitMassTauDown", &SVfitMassTauDown, &b_SVfitMassTauDown);
   fChain->SetBranchAddress("MCSignalParticle_p4", &MCSignalParticle_p4, &b_MCSignalParticle_p4);
   fChain->SetBranchAddress("MCSignalParticle_pdgid", &MCSignalParticle_pdgid, &b_MCSignalParticle_pdgid);
   fChain->SetBranchAddress("MCSignalParticle_charge", &MCSignalParticle_charge, &b_MCSignalParticle_charge);
   fChain->SetBranchAddress("MCSignalParticle_Poca", &MCSignalParticle_Poca, &b_MCSignalParticle_Poca);
   fChain->SetBranchAddress("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx, &b_MCSignalParticle_Tauidx);
   fChain->SetBranchAddress("MCTauandProd_p4", &MCTauandProd_p4, &b_MCTauandProd_p4);
   fChain->SetBranchAddress("MCTauandProd_Vertex", &MCTauandProd_Vertex, &b_MCTauandProd_Vertex);
   fChain->SetBranchAddress("MCTauandProd_pdgid", &MCTauandProd_pdgid, &b_MCTauandProd_pdgid);
   fChain->SetBranchAddress("MCTauandProd_midx", &MCTauandProd_midx, &b_MCTauandProd_midx);
   fChain->SetBranchAddress("MCTauandProd_charge", &MCTauandProd_charge, &b_MCTauandProd_charge);
   fChain->SetBranchAddress("MCTau_JAK", &MCTau_JAK, &b_MCTau_JAK);
   fChain->SetBranchAddress("MCTau_DecayBitMask", &MCTau_DecayBitMask, &b_MCTau_DecayBitMask);
   fChain->SetBranchAddress("MC_p4", &MC_p4, &b_MC_p4);
   fChain->SetBranchAddress("MC_pdgid", &MC_pdgid, &b_MC_pdgid);
   fChain->SetBranchAddress("MC_charge", &MC_charge, &b_MC_charge);
   fChain->SetBranchAddress("MC_midx", &MC_midx, &b_MC_midx);
   fChain->SetBranchAddress("MC_childpdgid", &MC_childpdgid, &b_MC_childpdgid);
   fChain->SetBranchAddress("MC_childidx", &MC_childidx, &b_MC_childidx);
   fChain->SetBranchAddress("MC_status", &MC_status, &b_MC_status);
   fChain->SetBranchAddress("SVfitMass", &SVfitMass, &b_SVfitMass);
   fChain->SetBranchAddress("SVfitTransverseMass", &SVfitTransverseMass, &b_SVfitTransverseMass);
   fChain->SetBranchAddress("SVfit_pt", &SVfit_pt, &b_SVfit_pt);
   fChain->SetBranchAddress("SVfit_ptUnc", &SVfit_ptUnc, &b_SVfit_ptUnc);
   fChain->SetBranchAddress("SVfit_eta", &SVfit_eta, &b_SVfit_eta);
   fChain->SetBranchAddress("SVfit_etaUnc", &SVfit_etaUnc, &b_SVfit_etaUnc);
   fChain->SetBranchAddress("SVfit_phi", &SVfit_phi, &b_SVfit_phi);
   fChain->SetBranchAddress("SVfit_phiUnc", &SVfit_phiUnc, &b_SVfit_phiUnc);
   fChain->SetBranchAddress("SVfit_fitMETRho", &SVfit_fitMETRho, &b_SVfit_fitMETRho);
   fChain->SetBranchAddress("SVfit_fitMETPhi", &SVfit_fitMETPhi, &b_SVfit_fitMETPhi);
   fChain->SetBranchAddress("isOSCand", &isOSCand, &b_isOSCand);
   fChain->SetBranchAddress("METx", &METx, &b_METx);
   fChain->SetBranchAddress("METy", &METy, &b_METy);
   fChain->SetBranchAddress("uncorrMETx", &uncorrMETx, &b_uncorrMETx);
   fChain->SetBranchAddress("uncorrMETy", &uncorrMETy, &b_uncorrMETy);
   fChain->SetBranchAddress("MET_cov00", &MET_cov00, &b_MET_cov00);
   fChain->SetBranchAddress("MET_cov01", &MET_cov01, &b_MET_cov01);
   fChain->SetBranchAddress("MET_cov10", &MET_cov10, &b_MET_cov10);
   fChain->SetBranchAddress("MET_cov11", &MET_cov11, &b_MET_cov11);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("mT_Dau1", &mT_Dau1, &b_mT_Dau1);
   fChain->SetBranchAddress("mT_Dau2", &mT_Dau2, &b_mT_Dau2);
   fChain->SetBranchAddress("PDGIdDaughters", &PDGIdDaughters, &b_PDGIdDaughters);
   fChain->SetBranchAddress("indexDau1", &indexDau1, &b_indexDau1);
   fChain->SetBranchAddress("indexDau2", &indexDau2, &b_indexDau2);
   fChain->SetBranchAddress("particleType", &particleType, &b_particleType);
   fChain->SetBranchAddress("discriminator", &discriminator, &b_discriminator);
   fChain->SetBranchAddress("daughters_muonID", &daughters_muonID, &b_daughters_muonID);
   fChain->SetBranchAddress("daughters_typeOfMuon", &daughters_typeOfMuon, &b_daughters_typeOfMuon);
   fChain->SetBranchAddress("dxy", &dxy, &b_dxy);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("dxy_innerTrack", &dxy_innerTrack, &b_dxy_innerTrack);
   fChain->SetBranchAddress("dz_innerTrack", &dz_innerTrack, &b_dz_innerTrack);
   fChain->SetBranchAddress("daughters_rel_error_trackpt", &daughters_rel_error_trackpt, &b_daughters_rel_error_trackpt);
   fChain->SetBranchAddress("SIP", &SIP, &b_SIP);
   fChain->SetBranchAddress("daughters_iseleBDT", &daughters_iseleBDT, &b_daughters_iseleBDT);
   fChain->SetBranchAddress("daughters_iseleWP80", &daughters_iseleWP80, &b_daughters_iseleWP80);
   fChain->SetBranchAddress("daughters_iseleWP90", &daughters_iseleWP90, &b_daughters_iseleWP90);
   fChain->SetBranchAddress("daughters_eleMVAnt", &daughters_eleMVAnt, &b_daughters_eleMVAnt);
   fChain->SetBranchAddress("daughters_eleMVA_HZZ", &daughters_eleMVA_HZZ, &b_daughters_eleMVA_HZZ);
   fChain->SetBranchAddress("daughters_passConversionVeto", &daughters_passConversionVeto, &b_daughters_passConversionVeto);
   fChain->SetBranchAddress("daughters_eleMissingHits", &daughters_eleMissingHits, &b_daughters_eleMissingHits);
   fChain->SetBranchAddress("daughters_iseleChargeConsistent", &daughters_iseleChargeConsistent, &b_daughters_iseleChargeConsistent);
   fChain->SetBranchAddress("daughters_eleCUTID", &daughters_eleCUTID, &b_daughters_eleCUTID);
   fChain->SetBranchAddress("decayMode", &decayMode, &b_decayMode);
   fChain->SetBranchAddress("tauID", &tauID, &b_tauID);
   fChain->SetBranchAddress("combreliso", &combreliso, &b_combreliso);
   fChain->SetBranchAddress("combreliso03", &combreliso03, &b_combreliso03);
   fChain->SetBranchAddress("daughters_depositR03_tracker", &daughters_depositR03_tracker, &b_daughters_depositR03_tracker);
   fChain->SetBranchAddress("daughters_depositR03_ecal", &daughters_depositR03_ecal, &b_daughters_depositR03_ecal);
   fChain->SetBranchAddress("daughters_depositR03_hcal", &daughters_depositR03_hcal, &b_daughters_depositR03_hcal);
   fChain->SetBranchAddress("daughters_decayModeFindingOldDMs", &daughters_decayModeFindingOldDMs, &b_daughters_decayModeFindingOldDMs);
   fChain->SetBranchAddress("daughters_SCeta", &daughters_SCeta, &b_daughters_SCeta);
   fChain->SetBranchAddress("againstElectronMVA5category", &againstElectronMVA5category, &b_againstElectronMVA5category);
   fChain->SetBranchAddress("againstElectronMVA5raw", &againstElectronMVA5raw, &b_againstElectronMVA5raw);
   fChain->SetBranchAddress("byPileupWeightedIsolationRaw3Hits", &byPileupWeightedIsolationRaw3Hits, &b_byPileupWeightedIsolationRaw3Hits);
   fChain->SetBranchAddress("footprintCorrection", &footprintCorrection, &b_footprintCorrection);
   fChain->SetBranchAddress("neutralIsoPtSumWeight", &neutralIsoPtSumWeight, &b_neutralIsoPtSumWeight);
   fChain->SetBranchAddress("photonPtSumOutsideSignalCone", &photonPtSumOutsideSignalCone, &b_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("daughters_decayModeFindingNewDMs", &daughters_decayModeFindingNewDMs, &b_daughters_decayModeFindingNewDMs);
   fChain->SetBranchAddress("daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits", &daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("daughters_byIsolationMVA3oldDMwoLTraw", &daughters_byIsolationMVA3oldDMwoLTraw, &b_daughters_byIsolationMVA3oldDMwoLTraw);
   fChain->SetBranchAddress("daughters_byIsolationMVA3oldDMwLTraw", &daughters_byIsolationMVA3oldDMwLTraw, &b_daughters_byIsolationMVA3oldDMwLTraw);
   fChain->SetBranchAddress("daughters_byIsolationMVA3newDMwoLTraw", &daughters_byIsolationMVA3newDMwoLTraw, &b_daughters_byIsolationMVA3newDMwoLTraw);
   fChain->SetBranchAddress("daughters_byIsolationMVA3newDMwLTraw", &daughters_byIsolationMVA3newDMwLTraw, &b_daughters_byIsolationMVA3newDMwLTraw);
   fChain->SetBranchAddress("daughters_byIsolationMVArun2v1DBoldDMwLTraw", &daughters_byIsolationMVArun2v1DBoldDMwLTraw, &b_daughters_byIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("daughters_chargedIsoPtSum", &daughters_chargedIsoPtSum, &b_daughters_chargedIsoPtSum);
   fChain->SetBranchAddress("daughters_neutralIsoPtSum", &daughters_neutralIsoPtSum, &b_daughters_neutralIsoPtSum);
   fChain->SetBranchAddress("daughters_puCorrPtSum", &daughters_puCorrPtSum, &b_daughters_puCorrPtSum);
   fChain->SetBranchAddress("daughters_numChargedParticlesSignalCone", &daughters_numChargedParticlesSignalCone, &b_daughters_numChargedParticlesSignalCone);
   fChain->SetBranchAddress("daughters_numNeutralHadronsSignalCone", &daughters_numNeutralHadronsSignalCone, &b_daughters_numNeutralHadronsSignalCone);
   fChain->SetBranchAddress("daughters_numPhotonsSignalCone", &daughters_numPhotonsSignalCone, &b_daughters_numPhotonsSignalCone);
   fChain->SetBranchAddress("daughters_daughters_numParticlesSignalCone", &daughters_daughters_numParticlesSignalCone, &b_daughters_daughters_numParticlesSignalCone);
   fChain->SetBranchAddress("daughters_numChargedParticlesIsoCone", &daughters_numChargedParticlesIsoCone, &b_daughters_numChargedParticlesIsoCone);
   fChain->SetBranchAddress("daughters_numNeutralHadronsIsoCone", &daughters_numNeutralHadronsIsoCone, &b_daughters_numNeutralHadronsIsoCone);
   fChain->SetBranchAddress("daughters_numPhotonsIsoCone", &daughters_numPhotonsIsoCone, &b_daughters_numPhotonsIsoCone);
   fChain->SetBranchAddress("daughters_numParticlesIsoCone", &daughters_numParticlesIsoCone, &b_daughters_numParticlesIsoCone);
   fChain->SetBranchAddress("daughters_leadChargedParticlePt", &daughters_leadChargedParticlePt, &b_daughters_leadChargedParticlePt);
   fChain->SetBranchAddress("daughters_trackRefPt", &daughters_trackRefPt, &b_daughters_trackRefPt);
   fChain->SetBranchAddress("daughters_isLastTriggerObjectforPath", &daughters_isLastTriggerObjectforPath, &b_daughters_isLastTriggerObjectforPath);
   fChain->SetBranchAddress("daughters_trgMatched", &daughters_trgMatched, &b_daughters_trgMatched);
   fChain->SetBranchAddress("daughters_isTriggerObjectforPath", &daughters_isTriggerObjectforPath, &b_daughters_isTriggerObjectforPath);
   fChain->SetBranchAddress("daughters_FilterFired", &daughters_FilterFired, &b_daughters_FilterFired);
   fChain->SetBranchAddress("daughters_isGoodTriggerType", &daughters_isGoodTriggerType, &b_daughters_isGoodTriggerType);
   fChain->SetBranchAddress("daughters_L3FilterFired", &daughters_L3FilterFired, &b_daughters_L3FilterFired);
   fChain->SetBranchAddress("daughters_L3FilterFiredLast", &daughters_L3FilterFiredLast, &b_daughters_L3FilterFiredLast);
   fChain->SetBranchAddress("daughters_HLTpt", &daughters_HLTpt, &b_daughters_HLTpt);
   fChain->SetBranchAddress("daughters_isL1IsoTau28Matched", &daughters_isL1IsoTau28Matched, &b_daughters_isL1IsoTau28Matched);
   fChain->SetBranchAddress("Muon_trackCharge", &Muon_trackCharge, &b_Muon_trackCharge);
   fChain->SetBranchAddress("Muon_pdgid", &Muon_pdgid, &b_Muon_pdgid);
   fChain->SetBranchAddress("Muon_B", &Muon_B, &b_Muon_B);
   fChain->SetBranchAddress("Muon_M", &Muon_M, &b_Muon_M);
   fChain->SetBranchAddress("Muon_par", &Muon_par, &b_Muon_par);
   fChain->SetBranchAddress("Muon_cov", &Muon_cov, &b_Muon_cov);
   fChain->SetBranchAddress("PFTauSVPos", &PFTauSVPos, &b_PFTauSVPos);
   fChain->SetBranchAddress("PFTauSVCov", &PFTauSVCov, &b_PFTauSVCov);
   fChain->SetBranchAddress("PFTauPionsP4", &PFTauPionsP4, &b_PFTauPionsP4);
   fChain->SetBranchAddress("PFTauPionsCharge", &PFTauPionsCharge, &b_PFTauPionsCharge);
   fChain->SetBranchAddress("PFTauSVChi2NDofMatchingQuality", &PFTauSVChi2NDofMatchingQuality, &b_PFTauSVChi2NDofMatchingQuality);
   fChain->SetBranchAddress("daughters_jetNDauChargedMVASel", &daughters_jetNDauChargedMVASel, &b_daughters_jetNDauChargedMVASel);
   fChain->SetBranchAddress("daughters_miniRelIsoCharged", &daughters_miniRelIsoCharged, &b_daughters_miniRelIsoCharged);
   fChain->SetBranchAddress("daughters_miniRelIsoNeutral", &daughters_miniRelIsoNeutral, &b_daughters_miniRelIsoNeutral);
   fChain->SetBranchAddress("daughters_jetPtRel", &daughters_jetPtRel, &b_daughters_jetPtRel);
   fChain->SetBranchAddress("daughters_jetPtRatio", &daughters_jetPtRatio, &b_daughters_jetPtRatio);
   fChain->SetBranchAddress("daughters_jetBTagCSV", &daughters_jetBTagCSV, &b_daughters_jetBTagCSV);
   fChain->SetBranchAddress("daughters_lepMVA_mvaId", &daughters_lepMVA_mvaId, &b_daughters_lepMVA_mvaId);
   fChain->SetBranchAddress("daughters_pca_x", &daughters_pca_x, &b_daughters_pca_x);
   fChain->SetBranchAddress("daughters_pca_y", &daughters_pca_y, &b_daughters_pca_y);
   fChain->SetBranchAddress("daughters_pca_z", &daughters_pca_z, &b_daughters_pca_z);
   fChain->SetBranchAddress("daughters_pcaRefitPV_x", &daughters_pcaRefitPV_x, &b_daughters_pcaRefitPV_x);
   fChain->SetBranchAddress("daughters_pcaRefitPV_y", &daughters_pcaRefitPV_y, &b_daughters_pcaRefitPV_y);
   fChain->SetBranchAddress("daughters_pcaRefitPV_z", &daughters_pcaRefitPV_z, &b_daughters_pcaRefitPV_z);
   fChain->SetBranchAddress("daughters_pcaGenPV_x", &daughters_pcaGenPV_x, &b_daughters_pcaGenPV_x);
   fChain->SetBranchAddress("daughters_pcaGenPV_y", &daughters_pcaGenPV_y, &b_daughters_pcaGenPV_y);
   fChain->SetBranchAddress("daughters_pcaGenPV_z", &daughters_pcaGenPV_z, &b_daughters_pcaGenPV_z);
   fChain->SetBranchAddress("JetsNumber", &JetsNumber, &b_JetsNumber);
   fChain->SetBranchAddress("PFTau_a1_lvp", &PFTau_a1_lvp, &b_PFTau_a1_lvp);
   fChain->SetBranchAddress("PFTau_a1_cov", &PFTau_a1_cov, &b_PFTau_a1_cov);
   fChain->SetBranchAddress("PFTau_a1_charge", &PFTau_a1_charge, &b_PFTau_a1_charge);
   fChain->SetBranchAddress("PFTau_a1_pdgid", &PFTau_a1_pdgid, &b_PFTau_a1_pdgid);
   fChain->SetBranchAddress("PFTau_a1_B", &PFTau_a1_B, &b_PFTau_a1_B);
   fChain->SetBranchAddress("PFTau_a1_M", &PFTau_a1_M, &b_PFTau_a1_M);

   fChain->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
   fChain->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
   fChain->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
   fChain->SetBranchAddress("jets_e", &jets_e, &b_jets_e);
   fChain->SetBranchAddress("jets_rawPt", &jets_rawPt, &b_jets_rawPt);
   fChain->SetBranchAddress("jets_area", &jets_area, &b_jets_area);
   fChain->SetBranchAddress("jets_mT", &jets_mT, &b_jets_mT);
   fChain->SetBranchAddress("jets_Flavour", &jets_Flavour, &b_jets_Flavour);
   fChain->SetBranchAddress("jets_HadronFlavour", &jets_HadronFlavour, &b_jets_HadronFlavour);
   fChain->SetBranchAddress("jets_genjetIndex", &jets_genjetIndex, &b_jets_genjetIndex);
   fChain->SetBranchAddress("jets_PUJetID", &jets_PUJetID, &b_jets_PUJetID);
   fChain->SetBranchAddress("jets_PUJetIDupdated", &jets_PUJetIDupdated, &b_jets_PUJetIDupdated);
   fChain->SetBranchAddress("jets_vtxPt", &jets_vtxPt, &b_jets_vtxPt);
   fChain->SetBranchAddress("jets_vtxMass", &jets_vtxMass, &b_jets_vtxMass);
   fChain->SetBranchAddress("jets_vtx3dL", &jets_vtx3dL, &b_jets_vtx3dL);
   fChain->SetBranchAddress("jets_vtxNtrk", &jets_vtxNtrk, &b_jets_vtxNtrk);
   fChain->SetBranchAddress("jets_vtx3deL", &jets_vtx3deL, &b_jets_vtx3deL);
   fChain->SetBranchAddress("jets_leadTrackPt", &jets_leadTrackPt, &b_jets_leadTrackPt);
   fChain->SetBranchAddress("jets_leptonPtRel", &jets_leptonPtRel, &b_jets_leptonPtRel);
   fChain->SetBranchAddress("jets_leptonPt", &jets_leptonPt, &b_jets_leptonPt);
   fChain->SetBranchAddress("jets_leptonDeltaR", &jets_leptonDeltaR, &b_jets_leptonDeltaR);
   fChain->SetBranchAddress("jets_chEmEF", &jets_chEmEF, &b_jets_chEmEF);
   fChain->SetBranchAddress("jets_chHEF", &jets_chHEF, &b_jets_chHEF);
   fChain->SetBranchAddress("jets_nEmEF", &jets_nEmEF, &b_jets_nEmEF);
   fChain->SetBranchAddress("jets_nHEF", &jets_nHEF, &b_jets_nHEF);
   fChain->SetBranchAddress("jets_MUF", &jets_MUF, &b_jets_MUF);
   fChain->SetBranchAddress("jets_neMult", &jets_neMult, &b_jets_neMult);
   fChain->SetBranchAddress("jets_chMult", &jets_chMult, &b_jets_chMult);
   fChain->SetBranchAddress("jets_jecUnc", &jets_jecUnc, &b_jets_jecUnc);
   fChain->SetBranchAddress("bDiscriminator", &bDiscriminator, &b_bDiscriminator);
   fChain->SetBranchAddress("bCSVscore", &bCSVscore, &b_bCSVscore);
   fChain->SetBranchAddress("pfCombinedMVAV2BJetTags", &pfCombinedMVAV2BJetTags, &b_pfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("PFjetID", &PFjetID, &b_PFjetID);
   fChain->SetBranchAddress("jetRawf", &jetRawf, &b_jetRawf);
   fChain->SetBranchAddress("ak8jets_px", &ak8jets_px, &b_ak8jets_px);
   fChain->SetBranchAddress("ak8jets_py", &ak8jets_py, &b_ak8jets_py);
   fChain->SetBranchAddress("ak8jets_pz", &ak8jets_pz, &b_ak8jets_pz);
   fChain->SetBranchAddress("ak8jets_e", &ak8jets_e, &b_ak8jets_e);
   fChain->SetBranchAddress("ak8jets_SoftDropMass", &ak8jets_SoftDropMass, &b_ak8jets_SoftDropMass);
   fChain->SetBranchAddress("ak8jets_PrunedMass", &ak8jets_PrunedMass, &b_ak8jets_PrunedMass);
   fChain->SetBranchAddress("ak8jets_TrimmedMass", &ak8jets_TrimmedMass, &b_ak8jets_TrimmedMass);
   fChain->SetBranchAddress("ak8jets_FilteredMass", &ak8jets_FilteredMass, &b_ak8jets_FilteredMass);
   fChain->SetBranchAddress("ak8jets_tau1", &ak8jets_tau1, &b_ak8jets_tau1);
   fChain->SetBranchAddress("ak8jets_tau2", &ak8jets_tau2, &b_ak8jets_tau2);
   fChain->SetBranchAddress("ak8jets_tau3", &ak8jets_tau3, &b_ak8jets_tau3);
   fChain->SetBranchAddress("ak8jets_CSV", &ak8jets_CSV, &b_ak8jets_CSV);
   fChain->SetBranchAddress("ak8jets_nsubjets", &ak8jets_nsubjets, &b_ak8jets_nsubjets);
   fChain->SetBranchAddress("subjets_px", &subjets_px, &b_subjets_px);
   fChain->SetBranchAddress("subjets_py", &subjets_py, &b_subjets_py);
   fChain->SetBranchAddress("subjets_pz", &subjets_pz, &b_subjets_pz);
   fChain->SetBranchAddress("subjets_e", &subjets_e, &b_subjets_e);
   fChain->SetBranchAddress("subjets_CSV", &subjets_CSV, &b_subjets_CSV);
   fChain->SetBranchAddress("subjets_ak8MotherIdx", &subjets_ak8MotherIdx, &b_subjets_ak8MotherIdx);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_cov", &pv_cov, &b_pv_cov);
   fChain->SetBranchAddress("pvRefit_x", &pvRefit_x, &b_pvRefit_x);
   fChain->SetBranchAddress("pvRefit_y", &pvRefit_y, &b_pvRefit_y);
   fChain->SetBranchAddress("pvRefit_z", &pvRefit_z, &b_pvRefit_z);
   fChain->SetBranchAddress("pvGen_x", &pvGen_x, &b_pvGen_x);
   fChain->SetBranchAddress("pvGen_y", &pvGen_y, &b_pvGen_y);
   fChain->SetBranchAddress("pvGen_z", &pvGen_z, &b_pvGen_z);
   fChain->SetBranchAddress("isRefitPV", &isRefitPV, &b_isRefitPV);
   Notify();
}

Bool_t NtupleReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupleReader_cxx
