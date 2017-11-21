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
