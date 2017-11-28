/*
 * DataMCCorrections.cxx
 *
 *  Created on: Nov  16, 2017
 *      Author: cherepanov
 */

#include "DataMCCorrections.h"
#include "Selection_Base.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

///////////////////////////
//
// Constructor
//

DataMCCorrections::DataMCCorrections(bool load_ZPtWeights){

	// define which scale factors should be loaded
	// individual SF can be switched off to avoid opening files which are not necessary
	loadZPtWeights = load_ZPtWeights;

	// define location of input root files
	TString basedir = "";
	basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

	// Z Pt Weights
	if(loadZPtWeights){
	  ZPtWeightFile= new TFile(basedir+"weights/zpt_weights_2016.root", "READ");
	  m_zPtHist = (TH2D*)(ZPtWeightFile->Get("zptmass_histo"));
	  m_zPtHistErr = (TH2D*)(ZPtWeightFile->Get("zptmass_histo_err"));
	}
	
} 

DataMCCorrections::~DataMCCorrections(){
 }

///////////////////////////
//
// Muon scale factors
//

float DataMCCorrections::ZPTWeight(float genMass, float genPt){
  float zPtReweight = m_zPtHist->GetBinContent(m_zPtHist->GetXaxis()->FindBin(genMass),m_zPtHist->GetYaxis()->FindBin(genPt));
  return zPtReweight;
}

float DataMCCorrections::ZPTWeightErr(float genMass, float genPt){
  float zPtReweightErr = m_zPtHistErr->GetBinContent(m_zPtHistErr->GetXaxis()->FindBin(genMass),m_zPtHistErr->GetYaxis()->FindBin(genPt));
  return zPtReweightErr;
}

float DataMCCorrections::AgainstElectronDataMCCorrection(TLorentzVector p4, TString type){
  float scale(1);
  if(type=="AgainstElectronMVAVLoose"){
    if(fabs(p4.Eta()) < 1.460)scale=1.213;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.375;
  }
  if(type=="AgainstElectronMVALoose"){
    if(fabs(p4.Eta()) < 1.460)scale=1.320;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.380;
  }
  if(type=="AgainstElectronMVAMedium"){
    if(fabs(p4.Eta()) < 1.460)scale=1.323;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.527;
  }
  if(type=="AgainstElectronMVATight"){
    if(fabs(p4.Eta()) < 1.460)scale=1.402;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.900;
  }
  if(type=="AgainstElectronMVAVTight"){
    if(fabs(p4.Eta()) < 1.460)scale=1.207;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.968;
  }
  return scale;
}


float DataMCCorrections::AgainstMuonDataMCCorrection(TLorentzVector p4, TString type){
  float scale(1);
  if(type=="AgainstMuonMVALoose3"){
    if(fabs(p4.Eta()) <0.4)scale=1.010;
    else if(fabs(p4.Eta()) > 0.4 && fabs(p4.Eta()) < 0.8) scale = 1.007;
    else if(fabs(p4.Eta()) > 0.8 && fabs(p4.Eta()) < 1.2) scale = 0.870;
    else if(fabs(p4.Eta()) > 1.2 && fabs(p4.Eta()) < 1.7) scale = 1.154;
    else if(fabs(p4.Eta()) > 1.7 && fabs(p4.Eta()) < 2.3) scale = 2.281;
  }
  if(type=="AgainstMuonMVATight3"){
    if(fabs(p4.Eta()) <0.4)scale=1.263;
    else if(fabs(p4.Eta()) > 0.4 && fabs(p4.Eta()) < 0.8) scale = 1.364;
    else if(fabs(p4.Eta()) > 0.8 && fabs(p4.Eta()) < 1.2) scale = 0.854;
    else if(fabs(p4.Eta()) > 1.2 && fabs(p4.Eta()) < 1.7) scale = 1.712;
    else if(fabs(p4.Eta()) > 1.7 && fabs(p4.Eta()) < 2.3) scale = 2.324;
  }
  return scale;
}

