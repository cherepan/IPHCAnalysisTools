#include "ZTauHTauH.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
 
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




ZTauHTauH::ZTauHTauH(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_pt(20),
  cMu_eta(2.1),
  cTau_pt(20),
  cTau_eta(2.1)
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ZTauHTauH::~ZTauHTauH(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTauHTauH::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)           cut.at(TriggerOk)=1;
    if(i==nGoodTaus)           cut.at(nGoodTaus)=2;
    if(i==FirstTauIsolation)   cut.at(FirstTauIsolation)=1;
    if(i==SecondTauIsolation)  cut.at(SecondTauIsolation)=1;
    if(i==nGoodMuons)          cut.at(nGoodMuons)=0;
    if(i==PairCharge)          cut.at(PairCharge)=0;
    if(i==PrimeVtx)            cut.at(PrimeVtx)=1;
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nGoodTaus){
      title.at(i)="Number of tau leptons ";
      hlabel="Number of tau leptons  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodTaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodTaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==FirstTauIsolation){
      title.at(i)="Isolation of First Tau ";
      hlabel="Isolation of First Tau  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_FirstTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_FirstTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==SecondTauIsolation){
      title.at(i)="Isolation of Second Tau  ";
      hlabel="Isolation of Second Tau  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SecondTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SecondTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nGoodMuons){
      title.at(i)="Number of muons ";
      hlabel="Number of muons  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==PairCharge){
      title.at(i)="Charge of a pair candidate";
      hlabel="Charge of a pair candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }



  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau candidate",30,24.5,80.5," P_{T}(#tau), GeV","Events");
  Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau candidate",30,24.5,99.5," E(#tau), GeV","Events");
  Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");

  Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau candidate",30,24.5,80.5," P_{T}(#tau), GeV","Events");
  Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau candidate",30,24.5,99.5," E(#tau), GeV","Events");
  Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");


  TauTauMass=HConfig.GetTH1D(Name+"_TauTauMass","Visible invariant mass of a tau pair",40,40,120," M(#tau#tau), GeV","Events");
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",6,0.5,6.5,"NQCD in ABCD","Events");

  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,0,2,"QCD Shape","");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",20,0,1," #Delta R","Events");

  Tau1Isolation=HConfig.GetTH1D(Name+"_Tau1Isolation","First Tau Isoaltion 1- Loose, 2- Medium, 3 Tight, 4-VTight",5,0.5,5.5," Discrimiantor","Events");
  Tau2Isolation=HConfig.GetTH1D(Name+"_Tau2Isolation","First Tau Isoaltion 1- Loose, 2- Medium, 3 Tight, 4-VTight",5,0.5,5.5," Discrimiantor","Events");


    Selection::ConfigureHistograms();   //   do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  ZTauHTauH::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&Tau1PT);
  Extradist1d.push_back(&Tau1E);
  Extradist1d.push_back(&Tau1HPSDecayMode);

  Extradist1d.push_back(&Tau2PT);
  Extradist1d.push_back(&Tau2E);
  Extradist1d.push_back(&Tau2HPSDecayMode);

  Extradist1d.push_back(&dRTauTau);
  Extradist1d.push_back(&TauTauMass);
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);
  Extradist1d.push_back(&Tau1Isolation);
  Extradist1d.push_back(&Tau2Isolation);
}

void  ZTauHTauH::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;

  
  //  int ntau(0); int nmu(0); 
  std::vector<int> goodMuonsIndex;
  std::vector<int> goodTauIndex;
  for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event
    if(Ntp->isTightGoodTau(iDaugther)){
      if(Ntp->tauBaselineSelection(iDaugther)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cTau_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cTau_eta  ){
	    goodTauIndex.push_back(iDaugther) ;  }}}}


    if(Ntp->isMediumGoodMuon(iDaugther)){
      if(Ntp->muonBaselineSelection(iDaugther)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cMu_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cMu_eta  ){
	    goodMuonsIndex.push_back(iDaugther) ;  }}}}
  }
    
 

  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;
  // std::cout<<"----------------   "<< std::endl;
  // std::cout<<" muon size   "<< goodMuonsIndex.size() <<std::endl;
  // for(unsigned int v = 0;v<goodMuonsIndex.size(); v++ ){
  //   std::cout<<"   muon index  "<<goodMuonsIndex.at(v)<< " pt   " <<Ntp->Daughters_P4(goodMuonsIndex.at(v)).Pt() << "  charge   "<<  Ntp->Daughters_charge(goodMuonsIndex.at(v))<< std::endl; 
  // }

  // std::cout<<" tau size   "<< goodTauIndex.size() <<std::endl;
  // for(unsigned int v = 0;v<goodTauIndex.size(); v++ ){
  //   std::cout<<"   tau index  "<<goodTauIndex.at(v)<< " pt   " <<Ntp->Daughters_P4(goodTauIndex.at(v)).Pt() << "  charge   "<<  Ntp->Daughters_charge(goodTauIndex.at(v))<< std::endl; 
  // }


  value.at(nGoodTaus)=goodTauIndex.size();
  pass.at(nGoodTaus) = (value.at(nGoodTaus) == cut.at(nGoodTaus));


  value.at(nGoodMuons)=goodMuonsIndex.size();
  pass.at(nGoodMuons) =(value.at(nGoodMuons) == cut.at(nGoodMuons));


  value.at(PairCharge) = ChargeSumDummy;
  value.at(FirstTauIsolation) = 0;
  value.at(SecondTauIsolation) = 0;
  if(goodTauIndex.size()==2){
    value.at(PairCharge) = Ntp->Daughters_charge(goodTauIndex.at(0)) + Ntp->Daughters_charge(goodTauIndex.at(1));
    value.at(FirstTauIsolation) = Ntp->isTightIsolatedTau(goodTauIndex.at(0));
    value.at(SecondTauIsolation) = Ntp->isTightIsolatedTau(goodTauIndex.at(1));

  }
  pass.at(PairCharge) = (value.at(PairCharge) == cut.at(PairCharge));
  pass.at(FirstTauIsolation) = (value.at(FirstTauIsolation) == cut.at(FirstTauIsolation));
  pass.at(SecondTauIsolation) = (value.at(SecondTauIsolation) == cut.at(SecondTauIsolation));

  // Here you can defined different type of weights you want to apply to events. At the moment only PU weight is considered if event is not data
  double wobs=1;
  double w;
  if(!Ntp->isData()){w = Ntp->PUReweight();}
  else{w=1;}

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
  exclude_cuts.push_back(FirstTauIsolation);
  exclude_cuts.push_back(SecondTauIsolation);
  exclude_cuts.push_back(PairCharge);
  std::cout<<" before  " << pass.at(TriggerOk) << "    " <<   pass.at(PrimeVtx) << "    " <<  pass.at(nGoodTaus)<< "    " <<   pass.at(FirstTauIsolation) << "    " <<  pass.at(SecondTauIsolation) << "    " <<  pass.at(nGoodMuons) << "    " <<  pass.at(PairCharge) << "  passAllBut  " << passAllBut(exclude_cuts) <<std::endl;

  if(passAllBut(exclude_cuts)){
    if(pass.at(FirstTauIsolation) && pass.at(SecondTauIsolation)){
      if(pass.at(PairCharge)){
	NQCD.at(t).Fill(1.,w); //A
      }  
      if(!pass.at(PairCharge)){
	NQCD.at(t).Fill(2.,w); //B
      }
      if(Ntp->isMediumIsolatedTau(goodTauIndex.at(0)) && Ntp->isLooseIsolatedTau(goodTauIndex.at(1))){
	if(pass.at(PairCharge)){
	  NQCD.at(t).Fill(3.,w); //С
	}
	if(!pass.at(PairCharge)){
	  NQCD.at(t).Fill(4.,w); //В
	}  
      }
    }
  }

  bool IsQCDEvent = false;
  if(!pass.at(PairCharge)){
    QCDShape.at(t).Fill(1,w);
    if(id == DataMCType::Data){
      t=HConfig.GetType(DataMCType::QCD);
      IsQCDEvent = true;
    }
  }
  if(IsQCDEvent)    pass.at(PairCharge)= true;
     

  std::vector<unsigned int> exclude_cuts_ForTauIso;
  exclude_cuts_ForTauIso.push_back(FirstTauIsolation);
  exclude_cuts_ForTauIso.push_back(SecondTauIsolation);

 if(passAllBut(exclude_cuts_ForTauIso)){
   if(Ntp->isLooseIsolatedTau(goodTauIndex.at(0)))Tau1Isolation.at(t).Fill(1.);
   if(Ntp->isMediumIsolatedTau(goodTauIndex.at(0)))Tau1Isolation.at(t).Fill(2.);
   if(Ntp->isTightIsolatedTau(goodTauIndex.at(0)))Tau1Isolation.at(t).Fill(3.);
   if(Ntp->isVTightIsolatedTau(goodTauIndex.at(0)))Tau1Isolation.at(t).Fill(4.);
   if(Ntp->isLooseIsolatedTau(goodTauIndex.at(1)))Tau2Isolation.at(t).Fill(1.);
   if(Ntp->isMediumIsolatedTau(goodTauIndex.at(1)))Tau2Isolation.at(t).Fill(2.);
   if(Ntp->isTightIsolatedTau(goodTauIndex.at(1)))Tau2Isolation.at(t).Fill(3.);
   if(Ntp->isVTightIsolatedTau(goodTauIndex.at(1)))Tau2Isolation.at(t).Fill(4.);


 }


  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){


 
  int TauIndex1 = goodTauIndex.at(0);
  int TauIndex2 = goodTauIndex.at(1);


  TLorentzVector Tau1P4 = Ntp->Daughters_P4(TauIndex1);
  TLorentzVector Tau2P4 = Ntp->Daughters_P4(TauIndex2);


  Tau1PT.at(t).Fill(Tau1P4.Pt(),w);  // Fill transverse momentum
  Tau1E.at(t).Fill(Tau1P4.E(),w);  // Fill transverse momentum
  Tau1HPSDecayMode.at(t).Fill(Ntp->decayMode(TauIndex1),w);

  Tau2PT.at(t).Fill(Tau2P4.Pt(),w);  // Fill transverse momentum
  Tau2E.at(t).Fill(Tau2P4.E(),w);  // Fill transverse momentum
  Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(TauIndex2),w);
  TauTauMass.at(t).Fill((Tau1P4+Tau2P4).M(),w);
  dRTauTau.at(t).Fill(Tau1P4.DeltaR(Tau2P4),w);
  }
}




//  This is a function if you want to do something after the event loop
void  ZTauHTauH::Finish(){
  if(mode == RECONSTRUCT){
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

    double OS2SS = (QCD_Integral_C_Data_minus_MC  - QCD_IntegralMC_C )/ (QCD_Integral_D_Data_minus_MC);
    double QCD_ScaleFactor = QCD_Integral_B_Data_minus_MC *OS2SS;


    std::cout << "OS/SS QCD Sample: " << OS2SS << std::endl;
    std::cout << "Scale Factor for QCD Sample: " << QCD_ScaleFactor << std::endl;
    std::cout << "QCD in B region "<<  QCD_Integral_B_Data_minus_MC <<std::endl;
    std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
    std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_ScaleFactor/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
  }

  Selection::Finish();
}






