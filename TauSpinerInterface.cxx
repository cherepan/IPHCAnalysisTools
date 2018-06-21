#include "TauSpinerInterface.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#ifdef USE_TauSpinner
#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TLorentzVector.h"

int TauSpinerInterface::signalcharge=-1;
bool TauSpinerInterface::initialized=false;
#endif

TauSpinerInterface::TauSpinerInterface(){
}

TauSpinerInterface::~TauSpinerInterface(){

}
#ifdef USE_TauSpinner
void TauSpinerInterface::Initialize(){
  Tauolapp::Tauola::initialize();
  Tauolapp::Tauola::setRandomGenerator(&TauSpinerInterface::CustomSeed);   // -----   change RANMAR seed
  //  std::cout<<"--------> seed is set to:    "<< TauSpinerInterface::CustomSeed() <<std::endl;
  string name="MSTW2008nnlo90cl.LHgrid";
  LHAPDF::initPDFSetByName(name);
  double CMSENE = 13000.0; // center of mass system energy.
                          // used in PDF calculation. For pp collisions only
  bool Ipp = true;  // for pp collisions 
  // Initialize TauSpinner
  //Ipol - polarization of input sample
  //nonSM2 - nonstandard model calculations
  //nonSMN
  int Ipol=2,nonSM2=0,nonSMN=0;

  initialize_spinner(Ipp,Ipol,nonSM2,nonSMN,CMSENE);

}



double TauSpinerInterface::Get(int type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2, ULong64_t EventNumber){
   if(!initialized){
     Initialize();
     initialized=true;  
   }
   //Initialize();

  double WT    = 1.0; // assume that there is 1 bosons decaying into                                                                                                                                                                        
  // Calculate weight for first boson
  if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 ){  
    WT = calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino
  }
  else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ){
    /*
    std::cout << "simpleParticle - 1 " << std::endl;
    std::cout << "simpleParticle " << X.pdgid() << " " << X.px() << " " <<X.py() << " " <<X.pz() << " " << X.e() << std::endl; 
    std::cout << "simpleParticle " <<  tau.pdgid() << " " << tau.px() << " " <<tau.py() << " " <<tau.pz() << " " << tau.e() 
	      << " size " << tau_daughters.size() << std::endl;
    TLorentzVector Tau(0,0,0,0);
    for(unsigned int i=0;i<tau_daughters.size();i++){
      std::cout << "simpleParticle" << tau_daughters.at(i).pdgid() << " " << tau_daughters.at(i).px() << " " 
		<<tau_daughters.at(i).py() << " " <<tau_daughters.at(i).pz() << " " << tau_daughters.at(i).e() << std::endl;
      Tau+=TLorentzVector(tau_daughters.at(i).px(),tau_daughters.at(i).py(),tau_daughters.at(i).pz(),tau_daughters.at(i).e());
    }
    std::cout << "Taureco " << Tau.Px() << " " << Tau.Py() << " " << Tau.Pz() << " " << Tau.E() << " " << Tau.M() <<  std::endl;
    std::cout << "simpleParticle - 2 " << std::endl;
    std::cout << "simpleParticle " << tau2.pdgid() << " " << tau2.px() << " " <<tau2.py() << " " <<tau2.pz() << " " << tau2.e() 
	      << " size " << tau_daughters2.size() <<std::endl;
    TLorentzVector Tau2(0,0,0,0);
    for(unsigned int i=0;i<tau_daughters2.size();i++){
      std::cout << "simpleParticle" << tau_daughters2.at(i).pdgid() << " " << tau_daughters2.at(i).px() << " "
		<<tau_daughters2.at(i).py() << " " <<tau_daughters2.at(i).pz() << " " << tau_daughters2.at(i).e() << std::endl;
      Tau2+=TLorentzVector(tau_daughters.at(i).px(),tau_daughters.at(i).py(),tau_daughters.at(i).pz(),tau_daughters.at(i).e());
    }
    std::cout << "Tau2reco " << Tau2.Px() << " " << Tau2.Py() << " " << Tau2.Pz() << " " << Tau2.E() << " " << Tau2.M() << std::endl;
    */



    int eventNumber = EventNumber;
    std::cout << "\nevent = " << EventNumber << std::endl;
    
    TauSpinner::SimpleParticle boson(9.057124, 4.738310, 362.649353, 373.200989, 23);
    
    TauSpinner::SimpleParticle tau1(-7.145843, 32.856289, 59.924850, 68.736778, -15);
    
    std::vector<TauSpinner::SimpleParticle> tauFinalStates1;
    tauFinalStates1.push_back(TauSpinner::SimpleParticle(-0.471688, 3.796974, 7.760280, 8.652251, -16)); // neutrino
    tauFinalStates1.push_back(TauSpinner::SimpleParticle(-6.674153, 29.059315, 52.164566, 60.084522, 211)); // pion
    
    TauSpinner::SimpleParticle tau2(16.202969, -28.117977, 302.724579, 304.464264, 15);
    
    std::vector<TauSpinner::SimpleParticle> tauFinalStates2;
    tauFinalStates2.push_back(TauSpinner::SimpleParticle(13.371717, -21.770473, 237.549759, 238.919739, 16)); // neutrino
    tauFinalStates2.push_back(TauSpinner::SimpleParticle(2.831248, -6.347505, 65.174820, 65.544518, -211)); // pion


    srand(391486);
    //    WT = calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
    WT = calculateWeightFromParticlesH(boson, tau1, tau2, tauFinalStates1,tauFinalStates2);
    double polSM=getTauSpin();
    std::cout<<"  TauSPinner Interface " <<std::endl;
    std::cout<<"  getTauSpin():   "<< polSM <<" WT:   "<<WT<<std::endl;
    // double WT    = 1.0;
    // double polSM=0.0;
    // int nIterations=100;
    // for (int iteration = 0; iteration < nIterations; ++iteration)
    //   {
    // 	WT = calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
	
    // 	// http://tauolapp.web.cern.ch/tauolapp/namespaceTauSpinner.html#af58cb8fff6c2c5bdd52f47c161fffe30
    // 	// http://tauolapp.web.cern.ch/tauolapp/tau__reweight__lib_8cxx_source.html#l00020
    // 	 polSM+= getTauSpin();
    //   }
    // polSM /= nIterations;




    //std::cout << "polSM=getTauSpin() " <<  polSM << " " << getTauSpin() << " WT " << WT << std::endl;
  if(type==hminus || type==hplus) WT=1.0;
    if(type==hminus && polSM>0.0) WT=0; // sign definition flipped in TauSpiner
    if(type==hplus && polSM<0.0) WT=0;
  }

  if(Spin==type || type==hplus ||  type==hminus) return WT;
  if(UnSpin==type) return 1.0/WT;
  if(FlipSpin==type) return (2.0-WT)/(WT);
  Logger(Logger::Warning) << "TauSpinerWeight TauSpinerType " << type << " is INVALID. Returning Spin WT=1.0." << std::endl;
  return 1.0;
}

double TauSpinerInterface::CustomSeed()
{
  double random;
  random = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  std::cout<<" random  seed: "<< random <<std::endl;
  return random;
}


#endif
