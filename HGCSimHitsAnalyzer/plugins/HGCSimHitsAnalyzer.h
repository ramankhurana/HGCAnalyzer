#ifndef HGCSimHitsAnalyzer_h
#define HGCSimHitsAnalyzer_h
// -*- C++ -*-
//
// Package:    HGCSimHitsAnalyzer
// Class:      HGCSimHitsAnalyzer
// 
/**\class HGCSimHitsAnalyzer HGCSimHitsAnalyzer.cc Validation/HGCalValidation/plugins/HGCSimHitsAnalyzer.cc

 Description: Validates SimHits of High Granularity Calorimeter

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raman Khurana
//         Created:  Fri, 31 Jan 2014 18:35:18 GMT
// $Id$
//
//


// system include files
#include <memory>
#include "TTree.h"
#include "../interface/hitsinfo.h"
#include "../interface/energysum.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include <CLHEP/Geometry/Transform3D.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include "TFile.h"
#include "TH3F.h"


//#include "LinkDef.h"

class HGCSimHitsAnalyzer : public edm::EDAnalyzer {

public:
  TTree* tree;

  unsigned int nlayers;
  float sampFracEE ;
  int nhits_[31];
  
  
  // SimHit Tree Vars
  //std::vector<std::pair<hitsinfo,energysum>>  SimHitVector;
  std::vector<hitsinfo>  SimHitVector;
  


  explicit HGCSimHitsAnalyzer(const edm::ParameterSet&);
  ~HGCSimHitsAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void analyzeHits (std::vector<PCaloHit>& hits);
  void FillHitsInfo(std::pair<hitsinfo,energysum> hit_, unsigned int itimeslice, double esum); 
  float OneOverMip(std::pair<hitsinfo,energysum> hit_);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  bool defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  std::string           nameDetector_, caloHitSource_, outputfilename_;
  DQMStore              *dbe_;
  HGCalDDDConstants     *hgcons_;
  int                   verbosity_;
  bool                  geometrydefined_, symmDet_;
  std::map<uint32_t, HepGeom::Transform3D> transMap_;
  TFile* f;
  TH1F* EnergyCalibrated_;
};


#endif
