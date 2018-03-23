#include "ZTauTau.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
//#include "SVfitProvider.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "SkimConfig.h"



#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TMath.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "Objects.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"



ZTauTau::ZTauTau(TString Name_, TString id_):
  Selection(Name_,id_),
  DataMC_Corr(true,true,false),
  tauTrgSF("tight")
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ZTauTau::~ZTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTauTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1;
    if(i==NPairsFound)         cut.at(NPairsFound)=1;
    if(i==Tau1Isolation)       cut.at(Tau1Isolation)=1.;
    if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    if(i==PairMass)            cut.at(PairMass)=100.;
    //if(i==MTM)                 cut.at(MTM)=40;
    
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // if(i==PrimeVtx){
    //   title.at(i)="Number of Prime Vertices $(N>$";
    //   title.at(i)+=cut.at(PrimeVtx);
    //   title.at(i)+=")";
    //   htitle=title.at(i);
    //   htitle.ReplaceAll("$","");
    //   htitle.ReplaceAll("\\","#");
    //   hlabel="Number of Prime Vertices";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    // }
    if(i==Trigger){
      title.at(i)="At least 1 good pair with Trig+Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NPairsFound){
      title.at(i)="Pairs with good DeltaR";
      hlabel="Pairs with good DeltaR";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     else if(i==Tau1Isolation){
      title.at(i)="Tau1 Isolation";
      hlabel="Isolation of Tau1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau2Isolation){
      title.at(i)="Tau2 Isolation";
      hlabel="Isolation of Tau2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton Veto  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,30,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,30,0,150,hlabel,"Events"));
      }

    /* else if(i==MTM){
      title.at(i)="Missing Transverse Mass";
      hlabel="Missing Transverse Mass";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTM_",htitle,30,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTM_",htitle,30,0,100,hlabel,"Events"));
      }*/
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau1 candidate",10,35,80," P_{T}(#tau1), GeV","Events");
  Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  Tau1Mass=HConfig.GetTH1D(Name+"_Tau1Mass","Mass of selected #tau1 candidate",10,0,2.," M(#tau1), GeV","Events");
  Tau1Phi=HConfig.GetTH1D(Name+"_Tau1Phi","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  Tau1Eta=HConfig.GetTH1D(Name+"_Tau1Eta","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  Tau1dz=HConfig.GetTH1D(Name+"_Tau1dz","Tau1 dz",10,-0.04,0.04,"Tau1 dz","Events");
  Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");

  Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau2 candidate",10,30,55," P_{T}(#tau2), GeV","Events");
  Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  Tau2Mass=HConfig.GetTH1D(Name+"_Tau2Mass","Mass of selected #tau2 candidate",10,0,2.," M(#tau2), GeV","Events");
  Tau2Phi=HConfig.GetTH1D(Name+"_Tau2Phi","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  Tau2Eta=HConfig.GetTH1D(Name+"_Tau2Eta","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  Tau2dz=HConfig.GetTH1D(Name+"_Tau2dz","Tau2dz",10,-0.04,0.04,"Tau2 dz","Events");
  Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");

  Tau1isolation=HConfig.GetTH1D(Name+"_Tau1isolation","First Tau isolation 1- Loose, 2- Medium, 3 Tight, 4-VTight",4,0.5,4.5," Discriminator","Events");
  Tau2isolation=HConfig.GetTH1D(Name+"_Tau2isolation","First Tau isolation 1- Loose, 2- Medium, 3 Tight, 4-VTight",4,0.5,4.5," Discriminator","Events");


  againstElectronVLooseMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronVLooseMVA6_Tau1","againstElectronVLooseMVA6_Tau1",2,-0.5,1.5,"againstElectronVLooseMVA6_Tau1","Events");
  againstElectronLooseMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronLooseMVA6_Tau1","againstElectronLooseMVA6_Tau1",2,-0.5,1.5,"againstElectronLooseMVA6_Tau1","Events");
  againstElectronMediumMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronMediumMVA6_Tau1","againstElectronMediumMVA6_Tau1",2,-0.5,1.5,"againstElectronMediumMVA6_Tau1","Events");
  againstElectronTightMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronTightMVA6_Tau1","againstElectronTightMVA6_Tau1",2,-0.5,1.5,"againstElectronTightMVA6_Tau1","Events");
  againstElectronVTightMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronVTightMVA6_Tau1","againstElectronVTightMVA6_Tau1",2,-0.5,1.5,"againstElectronVTightMVA6_Tau1","Events");
  againstMuonLoose3_Tau1=HConfig.GetTH1D(Name+"_againstMuonLoose3_Tau1","againstMuonLoose3_Tau1",2,-0.5,1.5,"againstMuonLoose3_Tau1","Events");
  againstMuonTight3_Tau1=HConfig.GetTH1D(Name+"_againstMuonTight3_Tau1","againstMuonTight3_Tau1",2,-0.5,1.5,"againstMuonTight3_Tau1","Events");
  byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1=HConfig.GetTH1D(Name+"_byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1","byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1",10,0,20,"byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1","Events");

  againstElectronVLooseMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronVLooseMVA6_Tau2","againstElectronVLooseMVA6_Tau2",2,-0.5,1.5,"againstElectronVLooseMVA6_Tau2","Events");
  againstElectronLooseMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronLooseMVA6_Tau2","againstElectronLooseMVA6_Tau2",2,-0.5,1.5,"againstElectronLooseMVA6_Tau2","Events");
  againstElectronMediumMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronMediumMVA6_Tau2","againstElectronMediumMVA6_Tau2",2,-0.5,1.5,"againstElectronMediumMVA6_Tau2","Events");
  againstElectronTightMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronTightMVA6_Tau2","againstElectronTightMVA6_Tau2",2,-0.5,1.5,"againstElectronTightMVA6_Tau2","Events");
  againstElectronVTightMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronVTightMVA6_Tau2","againstElectronVTightMVA6_Tau2",2,-0.5,1.5,"againstElectronVTightMVA6_Tau2","Events");
  againstMuonLoose3_Tau2=HConfig.GetTH1D(Name+"_againstMuonLoose3_Tau2","againstMuonLoose3_Tau2",2,-0.5,1.5,"againstMuonLoose3_Tau2","Events");
  againstMuonTight3_Tau2=HConfig.GetTH1D(Name+"_againstMuonTight3_Tau2","againstMuonTight3_Tau2",2,-0.5,1.5,"againstMuonTight3_Tau2","Events");
  byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2=HConfig.GetTH1D(Name+"_byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2","byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2",10,0,20,"byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2","Events");

  ExtraLeptonVeto=HConfig.GetTH1D(Name+"_ExtraLeptonVeto","ExtraLeptonVeto",2,-0.5,1.5,"ExtraLeptonVeto","Events");
  
  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,-0.5,1.5,"QCD Shape","");

  TauTauMass=HConfig.GetTH1D(Name+"_TauTauMass","Visible invariant mass of a tau pair",11,0,110," M(#tau#tau), GeV","Events");
  
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",7,0.,3.5," #Delta R","Events");

  MET=HConfig.GetTH1D(Name+"_MET","MET",20,0,80,"MET, GeV","Events");
  METphi=HConfig.GetTH1D(Name+"_METphi","METphi",10,-3.5,3.5,"METphi","Events");
  PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","PUPPImet",10,0,75,"PUPPImet, GeV","Events");
  PUPPImetphi=HConfig.GetTH1D(Name+"_PUPPImetphi","PUPPImetphi",10,-3.5,3.5,"PUPPImetphi","Events");
  TransverseMass=HConfig.GetTH1D(Name+"_TransverseMass","TransverseMass, GeV",22,0,110,"TransverseMass","Events");
  
  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,50,"N vtx","Events");
  NPU=HConfig.GetTH1D(Name+"_npu","npu",10,0,50,"N pu","Events");
  RHO=HConfig.GetTH1D(Name+"_rho","rho",10,0,30,"rho","Events");
  
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",12,0,12,"Number of jets","Events");

  h_SVFitMass = HConfig.GetTH1D(Name+"_SVFitMass","SVFitMass",100,0.,200.,"m_{SVfit}(#tau_{h},#tau_{h})/GeV");
  h_SVFitStatus = HConfig.GetTH1D(Name+"_SVFitStatus", "SVFitStatus", 5, -0.5, 4.5, "Status of SVFit calculation");

  svfTau1E = HConfig.GetTH1D(Name+"_svfTau1E","svFitTau1E",40,20.,120.,"E_{SVfit}(#tau_{h}1)/GeV");
  svfTau2E = HConfig.GetTH1D(Name+"_svfTau2E","svFitTau2E",40,20.,120.,"E_{SVfit}(#tau_{h}2)/GeV");

  Eta = HConfig.GetTH1D(Name+"_Eta", "Eta", 30, -5, 5, "Eta between Tau+ and Tau- in the new xy plan");

  Phi = HConfig.GetTH1D(Name+"_Phi","Phi",20,3.14,3.14,"Angle between Tau+ and an initial proton in the new xy plan");
  Theta = HConfig.GetTH1D(Name+"_Theta","Theta",20,3.14.,3.14,"Original theta of Tau-");
  
  /* CTN = HConfig.GetTH1D(Name+"_CTN","CTN",200,1.,-1.,"");
     CTT = HConfig.GetTH1D(Name+"_CTT","CTT",200,0.,2.,"");*/

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  ZTauTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&Tau1PT);
  Extradist1d.push_back(&Tau1E);
  Extradist1d.push_back(&Tau1Mass);
  Extradist1d.push_back(&Tau1Phi);
  Extradist1d.push_back(&Tau1Eta);
  Extradist1d.push_back(&Tau1dz);
  Extradist1d.push_back(&Tau1HPSDecayMode);
  
  Extradist1d.push_back(&Tau2PT);
  Extradist1d.push_back(&Tau2E);
  Extradist1d.push_back(&Tau2Mass);
  Extradist1d.push_back(&Tau2Phi);
  Extradist1d.push_back(&Tau2Eta);
  Extradist1d.push_back(&Tau2dz);
  Extradist1d.push_back(&Tau2HPSDecayMode);
  
  Extradist1d.push_back(&Tau1isolation);
  Extradist1d.push_back(&Tau2isolation);

  Extradist1d.push_back(&againstElectronVLooseMVA6_Tau1);
  Extradist1d.push_back(&againstElectronLooseMVA6_Tau1);
  Extradist1d.push_back(&againstElectronMediumMVA6_Tau1);
  Extradist1d.push_back(&againstElectronTightMVA6_Tau1);
  Extradist1d.push_back(&againstElectronVTightMVA6_Tau1);
  Extradist1d.push_back(&againstMuonLoose3_Tau1);
  Extradist1d.push_back(&againstMuonTight3_Tau1);
  Extradist1d.push_back(&byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1);

  Extradist1d.push_back(&againstElectronVLooseMVA6_Tau2);
  Extradist1d.push_back(&againstElectronLooseMVA6_Tau2);
  Extradist1d.push_back(&againstElectronMediumMVA6_Tau2);
  Extradist1d.push_back(&againstElectronTightMVA6_Tau2);
  Extradist1d.push_back(&againstElectronVTightMVA6_Tau2);
  Extradist1d.push_back(&againstMuonLoose3_Tau2);
  Extradist1d.push_back(&againstMuonTight3_Tau2);
  Extradist1d.push_back(&byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2);

  Extradist1d.push_back(&ExtraLeptonVeto);


  Extradist1d.push_back(&dRTauTau);
  Extradist1d.push_back(&TauTauMass);

  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);

  Extradist1d.push_back(&MET);
  Extradist1d.push_back(&METphi);
  Extradist1d.push_back(&PUPPImet);
  Extradist1d.push_back(&PUPPImetphi);
  Extradist1d.push_back(&TransverseMass);

  Extradist1d.push_back(&NPrimeVtx);
  Extradist1d.push_back(&NPU);
  Extradist1d.push_back(&RHO);


  Extradist1d.push_back(&NbJets);

Extradist1d.push_back(&h_SVFitMass);
 Extradist1d.push_back(&h_SVFitStatus);
 Extradist1d.push_back(&svfTau1E);
 Extradist1d.push_back(&svfTau2E);

  Extradist1d.push_back(&Eta);
  Extradist1d.push_back(&Phi);
  Extradist1d.push_back(&Theta);
}

void  ZTauTau::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;
  bool trig=0;
  std::vector<int> TauIndex ;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  value.at(Trigger)=0;
  MatchedTriggerNames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");
  MatchedTriggerNames.push_back("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
  /*for(int jj=0;jj<Ntp->NTriggers();jj++)
    {
    cout<<jj<<" "<<Ntp->TriggerName(jj)<<endl;
    }
    if(Ntp->TriggerAccept(TriggerIndexVector.at(0)))
    {
    if (Ntp->TriggerName(TriggerIndexVector.at(0)).Contains("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"))
    {
    trig0=1;
    }
    else if (Ntp->TriggerName(TriggerIndexVector.at(0)).Contains("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"))
    {
    trig1=1;
    }
    }
    if(Ntp->TriggerAccept(TriggerIndexVector.at(1)))
    {
      if (Ntp->TriggerName(TriggerIndexVector.at(1)).Contains("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"))
	{
	  trig0=1;
	}
      else if (Ntp->TriggerName(TriggerIndexVector.at(1)).Contains("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"))
	{
	  trig1=1;
	}
	}*/
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
      trig=1;
    }
  }
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
    if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),29) || Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),27))
      {
	TauIndex.push_back(iDaughter);
     }
  }
  if(trig && TauIndex.size()>1)value.at(Trigger)=1;
  pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger));
  value.at(Id_and_Kin)=0;
  int goodTau_counter=0;
  std::vector<int> thirdLeptonCounter;
  std::vector<int> goodTausIndex;
  for(unsigned int iDaughter=0;   iDaughter  < TauIndex.size()  ;iDaughter++ ) {
    if(Ntp->tauBaselineSelection(iDaughter,36., 2.1, 0,0)){
      goodTausIndex.push_back(TauIndex.at(iDaughter));
      goodTau_counter++;
      }
  }
  unsigned int Tau1= -1;
  unsigned int Tau2= -1;
  std::vector<int>  Sorted;
  std::vector<int>  PairsIndexTemp;
  std::vector<int>  PairsIndexTau1Temp;
  std::vector<int>  PairsIndexTau2Temp;
  int j=0;
  if(goodTau_counter>1)value.at(Id_and_Kin)=1;
  pass.at(Id_and_Kin)=(value.at(Id_and_Kin) == cut.at(Id_and_Kin));
  value.at(NPairsFound)=0;

  for(int unsigned ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(int unsigned jpair =1; jpair <goodTausIndex.size(); jpair++)
	{
	  if(jpair>ipair)
	    {
	      if((Ntp->Daughters_P4(goodTausIndex.at(ipair)).DeltaR(Ntp->Daughters_P4(goodTausIndex.at(jpair))))>0.5){
		PairsIndexTemp.push_back(j);
		PairsIndexTau1Temp.push_back(goodTausIndex.at(ipair));
		PairsIndexTau2Temp.push_back(goodTausIndex.at(jpair));
		j++;
	      }
	    }
	}
    }
  if(PairsIndexTemp.size()>0)value.at(NPairsFound)=1;
  pass.at(NPairsFound)=(value.at(NPairsFound)==cut.at(NPairsFound));
  if(pass.at(NPairsFound))
    {
      Sorted = Ntp->SortPair(PairsIndexTemp,PairsIndexTau1Temp,PairsIndexTau2Temp);
      Tau1=PairsIndexTau1Temp.at(Sorted.back());
      Tau2=PairsIndexTau2Temp.at(Sorted.back());
      value.at(Tau1Isolation)=0;
      value.at(Tau1Isolation) = (Ntp->isIsolatedTau(Tau1,"Tight"));
      pass.at(Tau1Isolation) = value.at(Tau1Isolation);
      value.at(Tau2Isolation)=0;
      value.at(Tau2Isolation) = (Ntp->isIsolatedTau(Tau2,"Tight"));
      pass.at(Tau2Isolation) = value.at(Tau2Isolation);

      
      value.at(LeptonVeto)=0;
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
	if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLeptonCounter.push_back(iDaughter);
	}
      }
      value.at(LeptonVeto) = thirdLeptonCounter.size()>0;
      pass.at(LeptonVeto) = (value.at(LeptonVeto)==cut.at(LeptonVeto));

      value.at(PairCharge)=0;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = 1;
      pass.at(PairCharge) = value.at(PairCharge);
      value.at(PairMass) = 999.;
      //value.at(MTM) = 999.;
      //value.at(MTM) = .;
      value.at(PairMass)=(Ntp->Daughters_P4(Tau1)+Ntp->Daughters_P4(Tau2)).M();
      pass.at(PairMass) = (value.at(PairMass) < cut.at(PairMass));
      //pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
    }
  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  if(!Ntp->isData() && id!=DataMCType::QCD) {
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
        //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
     if(!Ntp->isData() && pass.at(NPairsFound) ){
      double w1 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau1).Pt(),  Ntp->decayMode(Tau1)) ;  //from Luca
      double w2 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau2).Pt(),  Ntp->decayMode(Tau2)) ;
      w*=w1;
      w*=w2;
       }
      if(!Ntp->isData() && pass.at(NPairsFound) && id==33){
	w *= 0.95*0.95;
      }
  }
  TLorentzVector genMomentum(0,0,0,0);
  if( id == 33){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      // if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13)   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),5)  && Ntp->Genpart_status(imc) ==2){
      //   if(Ntp->Genpart_P4(imc).Pt() > 8){
      //     // std::cout<<" pdgid   "<< Ntp->Genpart_pdg(imc)<< "  px  " << Ntp->Genpart_P4(imc).Px() << " status   "<<Ntp->Genpart_status(imc) <<" index  "  << imc<<std::endl;
      //     // if(Ntp->Genpart_ZMothInd(imc)!=-1) 	std::cout<<"Mother    pdgid   "<< Ntp->Genpart_pdg(Ntp->Genpart_ZMothInd(imc))<< "  px  " << Ntp->Genpart_P4(Ntp->Genpart_ZMothInd(imc)).Px() <<std::endl;
      //   }
      // }
       if(fabs(Ntp->Genpart_pdg(imc)) ==15   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),0)&& Ntp->Genpart_status(imc) ==2) {
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  if( id == 30){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13) && Ntp->Genpart_status(imc) ==1  ){
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  
  float zptw(1);
  if(genMomentum.Pt()!=0 && genMomentum.M() > 75 && genMomentum.M() < 120){
    zptw = DataMC_Corr.ZPTWeight(genMomentum.M(),genMomentum.Pt());
  }
  w*=zptw;
  
  //-------------------------  mu/e tau fake rate weights 
  double wAgainstMuon1(1);
  double wAgainstElectron1(1);
  double wAgainstMuon2(1);
  double wAgainstElectron2(1);
  if(id == 33){
    if(pass.at(NPairsFound)){
      int matchedIndex1(-1);
      int matchedIndex2(-1);
      double dR1(999);
      double dR2(999);
      for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
	if(fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13)
	  {
	    if(sqrt(pow(Ntp->TauP4_Corrected(Tau1).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
		    pow(Ntp->TauP4_Corrected(Tau2).Eta() - Ntp->Genpart_P4(imc).Eta(),2)) < dR1)
	      {
		dR1 = sqrt(pow(Ntp->TauP4_Corrected(Tau1).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
			   pow(Ntp->TauP4_Corrected(Tau1).Eta() - Ntp->Genpart_P4(imc).Eta(),2));
		matchedIndex1=imc;
	      }
	    if(sqrt(pow(Ntp->TauP4_Corrected(Tau2).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
		    pow(Ntp->TauP4_Corrected(Tau2).Eta() - Ntp->Genpart_P4(imc).Eta(),2)) < dR2)
	      {
		dR2 = sqrt(pow(Ntp->TauP4_Corrected(Tau2).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
			   pow(Ntp->TauP4_Corrected(Tau2).Eta() - Ntp->Genpart_P4(imc).Eta(),2));
		matchedIndex2=imc;
	      }
	  }
      }
           if(dR1 < 0.2  && matchedIndex1!=-1 ){
	if(fabs(Ntp->Genpart_pdg(matchedIndex1)) ==13 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),5)   ) )
	  {
	    wAgainstMuon1 = DataMC_Corr.AgainstMuonDataMCCorrection(Ntp->TauP4_Corrected(Tau1),"AgainstMuonMVATight3");
	  }
	if(fabs(Ntp->Genpart_pdg(matchedIndex1)) ==11 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),5)   ) )
	  {
	    wAgainstElectron1 = DataMC_Corr.AgainstElectronDataMCCorrection(Ntp->TauP4_Corrected(Tau1),"AgainstElectronMVATight");
	  }
      }
           if(dR2 < 0.2  && matchedIndex2!=-1 ){
	if(fabs(Ntp->Genpart_pdg(matchedIndex2)) ==13 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),5)   ) )
	  {
	    wAgainstMuon2 = DataMC_Corr.AgainstMuonDataMCCorrection(Ntp->TauP4_Corrected(Tau2),"AgainstMuonMVATight3");
	  }
	if(fabs(Ntp->Genpart_pdg(matchedIndex2)) ==11 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),5)   ) )
	  {
	    wAgainstElectron2 = DataMC_Corr.AgainstElectronDataMCCorrection(Ntp->TauP4_Corrected(Tau2),"AgainstElectronMVATight");
	  }
      }
    }
  }
  w*=wAgainstMuon1;//cout<<"againstmu1: "<<w<<endl;
  w*=wAgainstElectron1;//cout<<"againstele1: "<<w<<endl;
  w*=wAgainstMuon2;//cout<<"againstmu2: "<<w<<endl;
  w*=wAgainstElectron2;//cout<<"againstele2: "<<w<<endl;
  //w*=Ntp->MC_weight();cout<<"MC"<<w<<endl; //generator weight

  // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  SS
   *        |        *       S
   * ----------------*------ i
   *        |        *       g
   *    A   |    B   *  OS   n
   *        |        *
   *******************
   *  Iso   | AntiIso
   *
   *     TauIsolation
   */

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(Tau1Isolation);
  exclude_cuts.push_back(Tau2Isolation);
  exclude_cuts.push_back(PairCharge);
  // std::cout<<" before  " << pass.at(TriggerOk) << "    " <<   pass.at(PrimeVtx) << "    " <<  pass.at(NPairsFound)<< "    " <<   pass.at(FirstTauIsolation) << "    " <<  pass.at(SecondTauIsolation) << "    " <<  pass.at(nGoodMuons) << "    " <<  pass.at(PairCharge) << "  passAllBut  " << passAllBut(exclude_cuts) <<std::endl;

  if(passAllBut(exclude_cuts)){
    //    for(unsigned int ia=0; ia<pass.size(); ia++){         std::cout<<" ia  "<< ia <<  "   pass  " <<pass.at(ia) << std::endl;    }
    // if(pass.at(FirstTauIsolation) && pass.at(SecondTauIsolation)){
    //   if(pass.at(PairCharge)){
    // 	NQCD.at(t).Fill(1.,w); //A
    //   }  
    //   if(!pass.at(PairCharge)){
    // 	NQCD.at(t).Fill(2.,w); //B
    //   }
    //   if(Ntp->isIsolatedTau(TauIndex_1,"Medium") && Ntp->isIsolatedTau(TauIndex_2,"Loose")){
    // 	if(pass.at(PairCharge)){
    // 	  NQCD.at(t).Fill(3.,w); //С
    // 	}
    // 	if(!pass.at(PairCharge)){
    // 	  NQCD.at(t).Fill(4.,w); //В
    // 	}
    //   }
    // }
    if(pass.at(PairCharge)){
      if(pass.at(Tau1Isolation) && pass.at(Tau2Isolation) ){
	NQCD.at(t).Fill(1.,w); //A
      }
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	NQCD.at(t).Fill(2.,w); //B
      }
    }
    if(!pass.at(PairCharge)){
      if(pass.at(Tau1Isolation) && pass.at(Tau2Isolation)){
	NQCD.at(t).Fill(3.,w); //C
      }
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	NQCD.at(t).Fill(4.,w); //D
      }
    }
  }
  bool IsQCDEvent = false;
  if(passAllBut(exclude_cuts)){
    if(pass.at(PairCharge)){
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	if(id == DataMCType::Data){
	  QCDShape.at(t).Fill(1,w);
	  t=HConfig.GetType(DataMCType::QCD);
	  IsQCDEvent = true;
	}
      }
    }
  }

  if(IsQCDEvent){ pass.at(PairCharge)= true;pass.at(Tau2Isolation)= true;pass.at(Tau1Isolation)=true;}
  
  std::vector<unsigned int> exclude_cuts_ForTauIso;
  exclude_cuts_ForTauIso.push_back(Tau1Isolation);
  exclude_cuts_ForTauIso.push_back(Tau2Isolation);
  if(passAllBut(exclude_cuts_ForTauIso)){
    if(Ntp->isIsolatedTau(Tau1,"Loose"))Tau1isolation.at(t).Fill(1.);
    if(Ntp->isIsolatedTau(Tau1,"Medium"))Tau1isolation.at(t).Fill(2.);
    if(Ntp->isIsolatedTau(Tau1,"Tight"))Tau1isolation.at(t).Fill(3.);
    if(Ntp->isIsolatedTau(Tau1,"VTight"))Tau1isolation.at(t).Fill(4.);
    if(Ntp->isIsolatedTau(Tau2,"Loose"))Tau2isolation.at(t).Fill(1.);
    if(Ntp->isIsolatedTau(Tau2,"Medium"))Tau2isolation.at(t).Fill(2.);
    if(Ntp->isIsolatedTau(Tau2,"Tight"))Tau2isolation.at(t).Fill(3.);
    if(Ntp->isIsolatedTau(Tau2,"VTight"))Tau2isolation.at(t).Fill(4.);
    }

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events which passed selection
  if(status) {
    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
    if(id !=DataMCType::Data && id !=DataMCType::QCD)	  pvx = Ntp->PUNumInteractions();
    NPrimeVtx.at(t).Fill(pvx,w);
    NPU.at(t).Fill(Ntp->npu(),w);
    RHO.at(t).Fill(Ntp->rho(),w);
  
    TLorentzVector Tau1P4 = Ntp->TauP4_Corrected(Tau1);
    TLorentzVector Tau2P4 = Ntp->TauP4_Corrected(Tau2);
    std::vector<int> thirdLepton;


    //---------  svfit ---------------------
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    classic_svFit::MeasuredTauLepton lep1(1, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M());
    classic_svFit::MeasuredTauLepton lep2(1, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M());

    measuredTauLeptons.push_back(lep1);
    measuredTauLeptons.push_back(lep2);
    TMatrixD metcov(2,2);
    double metx = Ntp->MET()*cos(Ntp->METphi());
    double mety = Ntp->MET()*sin(Ntp->METphi());

    metcov[0][0] = Ntp->PFMETCov00();
    metcov[1][0] = Ntp->PFMETCov01();
    metcov[0][1] = Ntp->PFMETCov10();
    metcov[1][1] = Ntp->PFMETCov11();

    svfitAlgo1.addLogM_fixed(true,5.0);
    svfitAlgo1.setDiTauMassConstraint(-1.0);
    svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );

    if(svfitAlgo1.isValidSolution()){
      double higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
      h_SVFitMass.at(t).Fill(higgsmass,w); 
    }
   
    ClassicSVfit svfitAlgo2;
    svfitAlgo2.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
    svfitAlgo2.addLogM_fixed(true, 5.);
    svfitAlgo2.integrate(measuredTauLeptons,metx,mety, metcov );

    classic_svFit::LorentzVector tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau1LV();
    classic_svFit::LorentzVector tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau2LV();

    svfTau1E.at(t).Fill(tau1P4.E(),w);
    svfTau2E.at(t).Fill(tau2P4.E(),w);
    
    Tau1PT.at(t).Fill(Tau1P4.Pt(),w);
    Tau1E.at(t).Fill(Tau1P4.E(),w);
    Tau1Mass.at(t).Fill(Tau1P4.M(),w);
    Tau1Phi.at(t).Fill(Tau1P4.Phi(),w);
    Tau1Eta.at(t).Fill(Tau1P4.Eta(),w);
    Tau1dz.at(t).Fill(Ntp->dz(Tau1),w);
    Tau1HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau1),w);

    Tau2PT.at(t).Fill(Tau2P4.Pt(),w);
    Tau2E.at(t).Fill(Tau2P4.E(),w);
    Tau2Mass.at(t).Fill(Tau2P4.M(),w);
    Tau2Phi.at(t).Fill(Tau2P4.Phi(),w);
    Tau2Eta.at(t).Fill(Tau2P4.Eta(),w);
    Tau2dz.at(t).Fill(Ntp->dz(Tau2),w);
    Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau2),w);
    againstElectronVLooseMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronVLooseMVA6),w);
    againstElectronLooseMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronLooseMVA6),w);
    againstElectronMediumMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronMediumMVA6),w);
    againstElectronTightMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronTightMVA6),w);
    againstElectronVTightMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronVTightMVA6),w);
    againstMuonLoose3_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstMuonLoose3),w);
    againstMuonTight3_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstMuonTight3),w);
    byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1.at(t).Fill(Ntp->Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(Tau1),w);

    againstElectronVLooseMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronVLooseMVA6),w);
    againstElectronLooseMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronLooseMVA6),w);
    againstElectronMediumMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronMediumMVA6),w);
    againstElectronTightMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronTightMVA6),w);
    againstElectronVTightMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronVTightMVA6),w);
    againstMuonLoose3_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstMuonLoose3),w);
    againstMuonTight3_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstMuonTight3),w);
    byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2.at(t).Fill(Ntp->Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(Tau2),w);
    
    for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
      if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLepton.push_back(iDaughter);
      }
    }
    if(thirdLepton.size()>0)ExtraLeptonVeto.at(t).Fill(1.,w);
    else ExtraLeptonVeto.at(t).Fill(0.,w);
      
    TauTauMass.at(t).Fill((Tau1P4+Tau2P4).M(),w);
    dRTauTau.at(t).Fill(Tau1P4.DeltaR(Tau2P4),w);

    MET.at(t).Fill(Ntp->MET(),w);
    METphi.at(t).Fill(Ntp->METphi(),w);
    PUPPImet.at(t).Fill(Ntp->PUPPImet(),w);
    PUPPImetphi.at(t).Fill(Ntp->PUPPImetphi(),w);
    TransverseMass.at(t).Fill(Ntp->transverseMass(Tau1P4.Pt(), Tau1P4.Phi(), Tau2P4.Pt(), Tau2P4.Phi()),w);

  
    int jets_counter=0;
    for(int ijet=0; ijet< Ntp->JetsNumber(); ijet++) {
      if((((Ntp->Jet_P4(ijet)).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<4.7))) {
	if((((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Tau1)))>0.5) && (((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Tau2)))>0.5))jets_counter++; {
	  //if(((Ntp->Jet_P4(ijet).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<2.4)) && (((Ntp->Jet_P4(ijet))->(Ntp->bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags'))) > 0.800)) {
	  if((((Ntp->jets_nHEF(ijet))<0.99) && ((Ntp->jets_nEmEF(ijet)<0.99)) && ((Ntp->NumConst(ijet))>1)) && (((fabs((Ntp->Jet_P4(ijet).Eta())<=2.4)) && ((Ntp->jets_chHEF(ijet))>0) && ((Ntp->jets_chMult(ijet))>0) && ((Ntp->jets_chEmEF(ijet))<0.99)) || (fabs((Ntp->Jet_P4(ijet).Eta())>2.4))) && (fabs((Ntp->Jet_P4(ijet)).Eta()<=2.7)))jets_counter++;
	  else if(((Ntp->jets_nHEF(ijet))<0.98) && ((Ntp->jets_nEmEF(ijet))>0.01) && ((Ntp->jets_neMult(ijet))>2) && (fabs((Ntp->Jet_P4(ijet).Eta())>2.7)) && (fabs((Ntp->Jet_P4(ijet).Eta())<=3.0)))jets_counter++;
	  else if(((Ntp->jets_nEmEF(ijet))<0.90) && ((Ntp->jets_neMult(ijet))>10) && (fabs((Ntp->Jet_P4(ijet).Eta())>3.0)))jets_counter++;
	}
      }
    }

    NbJets.at(t).Fill(jets_counter,w);

    classic_svFit::LorentzVector Tauplus;
    classic_svFit::LorentzVector Tauminus;
    if(Ntp->Daughters_charge(Tau1)>0)
      {
	Tauplus=tau1P4;
	Tauminus=tau2P4;
      }
    else
      {
	Tauplus=tau2P4;
	Tauminus=tau1P4;
      }
    
    TVector3 Z(Tauminus);
    TVector3 Y=Tauminus.Cross(Tauplus);
    TVector3 X=Y.Cross(Tauminus);
    TVector3 tauplus=(Tauplus*X,Tauplus*Y,Tauplus*Z);
    TVector3 tauminus=(Tauminus*X,Tauminus*Y,Tauminus*Z);
    TVector3 proton=(0,0,1);
    TVector3 xyproton=(proton*X,proton*Y,0);
    TVector3 xytauplus=(tauplus*X,tauplus*Y,0);
    
    Eta.at(t).Fill(TMath::ATanH((tauplus.Z())/(sqrt((tauplus.Z())*(tauplus.Z())+(tauplus.Y())*(tauplus.Y())+(tauplus.X())*(tauplus.X())))),w);
    Phi.at(t).Fill(xytauplus.Angle(xyproton),w);
    Theta.at(t).Fill(Tauminus.Theta(),w);
  }
}
//  This is a function if you want to do something after the event loop
void  ZTauTau::Finish() {
  if(mode == RECONSTRUCT) {
      std::cout<<" Starting Finish!  " <<std::endl;
    
      std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
      std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
      std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
      std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
      SkimConfig SC;
      SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);

      std::vector<double> QCD_Integral_B;
      double QCD_IntegralMC_B;
      double QCD_Integral_B_Data_minus_MC = 0;
    
      std::vector<double> QCD_Integral_C;
      double QCD_IntegralMC_C;
      double QCD_Integral_C_Data_minus_MC = 0;
    
      std::vector<double> QCD_Integral_D;
      double QCD_IntegralMC_D;
      double QCD_Integral_D_Data_minus_MC = 0;

      //Get Yields in ABCD for QCD Scalefactor                                                                                                                                                                  
      for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
      QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
      QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
      if(CrossSectionandAcceptance.at(i)>0){
      QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      }
      }
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(HConfig.GetID(i) == DataMCType::Data){
	QCD_Integral_B_Data_minus_MC  += QCD_Integral_B.at(i);
	QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
	QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
      }
      if(CrossSectionandAcceptance.at(i)>0){
	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
      }
    }

    double CDFactor = (QCD_Integral_C_Data_minus_MC  - QCD_IntegralMC_C )/ (QCD_Integral_D_Data_minus_MC - QCD_IntegralMC_D);
    double QCD_Signal = QCD_Integral_B_Data_minus_MC *CDFactor;


    std::cout << "Factor: " << CDFactor << std::endl;
    std::cout << "QCD_Signal: " << QCD_Signal << std::endl;
    std::cout << "QCD in B region "<<  QCD_Integral_B_Data_minus_MC <<std::endl;
    std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
    std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_Signal/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    }
  Selection::Finish();
}
