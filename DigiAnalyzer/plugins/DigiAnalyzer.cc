// -*- C++ -*-
//
// Package:    DigiAnalyzer
// Class:      DigiAnalyzer
// 
/**\class DigiAnalyzer DigiAnalyzer.cc HGCAnalyzer/DigiAnalyzer/plugins/DigiAnalyzer.cc

This class saves the digi information in a tree for further usage. 
Digi information include : 
    ADC counts
    Gain 
    Charge 
    

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raman Khurana
//         Created:  Thu, 04 Sep 2014 15:37:46 GMT
// $Id$
//
//


// system include files
#include <memory>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Validation/HGCalValidation/plugins/HGCalDigiValidation.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
/////                                                                                                                                                                  
#include "DetectorDescription/Core/interface/DDExpandedView.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "FWCore/Utilities/interface/InputTag.h" 
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include <cmath>

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#define SAMPLES 3

class DigiAnalyzer : public edm::EDAnalyzer {
public:
  
  edm::Service<TFileService> fs;
  TTree* digiTree;
  TH1F* ADC_[3];
  std::vector<int> ADC_All;
  std::vector<int> ADC_TS1;
  std::vector<int> ADC_TS2;
  std::vector<int> ADC_TS3;
  std::vector<int> Gain;
  std::string  nameDetector_, DigiSource_;
  HGCalDDDConstants     *hgcons_;
  explicit DigiAnalyzer(const edm::ParameterSet&);
  ~DigiAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  template<class T1, class T2> 
  void HGCDigiSaver(T1 detId, ForwardSubdetector subdet, const HGCalGeometry& geom0, const T2 it);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DigiAnalyzer::DigiAnalyzer(const edm::ParameterSet& iConfig)
{
  nameDetector_  = iConfig.getParameter<std::string>("DetectorName");
  DigiSource_ = iConfig.getParameter<std::string>("DigiSource");
  
  fs->make<TTree>("tree","tree");
  digiTree = fs->make<TTree>("digiTree","digiTree");
  
  for(int i=0; i<SAMPLES; i++){
    TString postname = Form("_%d",i);
    ADC_[i] = fs->make<TH1F> ("ADC"+postname, "ADC", 100,0,100);
  }
  //now do what ever initialization is needed

}


DigiAnalyzer::~DigiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  ADC_All.clear();
  ADC_TS1.clear();
  ADC_TS2.clear();
  ADC_TS3.clear();
  Gain.clear();
  
  edm::ESHandle<HGCalGeometry> geom;
  ForwardSubdetector subdet;
  
  if (nameDetector_ == "HGCalEESensitive") {
    subdet = HGCEE;
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geom);
    
    if (!geom.isValid())
      std::cout << "Cannot get valid HGCalGeometry Object for " << "HGCalEESensitive" << std::endl;
    
    const HGCalGeometry& geom0 = *geom;
    edm::Handle<HGCEEDigiCollection> theHGCEEDigiContainers;
    iEvent.getByLabel("mix", DigiSource_, theHGCEEDigiContainers);
    
    for(HGCEEDigiCollection::const_iterator it =theHGCEEDigiContainers->begin();
	it !=theHGCEEDigiContainers->end(); ++it) {
      HGCEEDetId detId = it->id();
      HGCDigiSaver(detId, subdet, geom0, it);
    }
  }

  else if((nameDetector_ == "HGCalHESiliconSensitive") || (nameDetector_ == "HGCalHEScintillatorSensitive")){
    
    if (nameDetector_ == "HGCalHESiliconSensitive") subdet = HGCHEF;
    else subdet = HGCHEB;
    
    iSetup.get<IdealGeometryRecord>().get(nameDetector_,geom);
    edm::Handle<HGCHEDigiCollection> theHGCHEDigiContainers;
    iEvent.getByLabel("mix", DigiSource_, theHGCHEDigiContainers);
    
    for(HGCHEDigiCollection::const_iterator it =theHGCHEDigiContainers->begin();
	it !=theHGCHEDigiContainers->end(); ++it) {
      HGCHEDetId detId = it->id();
      
      if (!geom.isValid())
	std::cout << "Cannot get valid HGCalGeometry Object for " << "HGCalEESensitive" << std::endl;
      
      const HGCalGeometry& geom0 = *geom;
      HGCDigiSaver(detId, subdet, geom0, it);
    }
  }
  else
    std::cout << "invalid detector name !!" << std::endl;

  digiTree->Fill();
}



template<class T1, class T2> 
void DigiAnalyzer::HGCDigiSaver(T1 detId, ForwardSubdetector subdet, const HGCalGeometry& geom0, const T2 it) {
  double digiCharge = 0;
  int nSample = it->size();
  for (int i=0; i<nSample; ++i) {
    HGCSample hgcSample = it->sample(i);
    uint16_t gain = hgcSample.gain();
    uint16_t adc = hgcSample.adc();
    double charge = adc*gain;
    ADC_[i]->Fill((int)adc);
    ADC_All.push_back((int)adc);
    if(i==0) ADC_TS1.push_back((int)adc);
    if(i==1) ADC_TS2.push_back((int)adc);
    if(i==2) ADC_TS3.push_back((int)adc);
    Gain.push_back((int)gain);
    digiCharge += charge;
  }
  
}



// ------------ method called once each job just before starting event loop  ------------
void 
DigiAnalyzer::beginJob()
{
  digiTree->Branch("ADC_All", "std::vector<int>", &ADC_All);
  digiTree->Branch("ADC_TS1", "std::vector<int>", &ADC_TS1);
  digiTree->Branch("ADC_TS2", "std::vector<int>", &ADC_TS2);
  digiTree->Branch("ADC_TS3", "std::vector<int>", &ADC_TS3);
  digiTree->Branch("Gain","std::vector<int>",&Gain);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DigiAnalyzer::endJob() 
{
  
  //digiTree->Write();
  //for(int i=0; i<SAMPLES; i++){
  //  ADC_[i]->Write();
  //}
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DigiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DigiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DigiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DigiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiAnalyzer);
