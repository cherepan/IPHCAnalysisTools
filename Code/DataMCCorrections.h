/*
 * DataMCCorrections.cxx
 *
 *  Created on: Nov  16, 2017
 *      Author: cherepanov
 */


#ifndef DataMCCorrections_H_
#define DataMCCorrections_H_

#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"

class DataMCCorrections {

public:

	DataMCCorrections( bool load_ZPtWeights = true);
	virtual ~DataMCCorrections();


	// Higgs/Z pT reweighting
	double HiggsPtWeight(TLorentzVector vect, int mass, TString shift = "nominal");
	float ZPTWeight(float genMass, float genPt);
	float ZPTWeightErr(float genMass, float genPt);

	float AgainstElectronDataMCCorrection(TLorentzVector p4, TString type);
	float AgainstMuonDataMCCorrection(TLorentzVector p4, TString type);
private:
	// booleans to switch on/off individual scale factors
	bool loadZPtWeights;
	//
	// Root files for scale factors
	//

	// Z pT reweighting
	TFile* ZPtWeightFile;

	//
	// Histograms for scale factors
	//

	// Z pT reweighting
	TH2D* m_zPtHist;
	TH2D* m_zPtHistErr;

};



#endif /* REFERENCESCALEFACTORS_H_ */
