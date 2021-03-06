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
//#include "Validation/HGCalValidation/plugins/HGCalDigiValidation.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include <CLHEP/Geometry/Transform3D.h>


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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <algorithm>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "../interface/IsolationComputation.h"
#define SAMPLES 6
#define nlayer 31
#define ncell 5
class DigiAnalyzer : public edm::EDAnalyzer {
public:
  
  edm::Service<TFileService> fs;
  TTree* digiTree;
  HGCalTopology* hgctopology_;
  IsolationComputation* isocomputation;
  ForwardSubdetector subdet;
  typedef std::map<HGCEEDetId, int> ADCInLayer;
  ADCInLayer AllADCInLayer; ;
  std::array<ADCInLayer,31> AllADCInDet;
  int maximumCell_;
  TH1F* SeedCell_[nlayer];
  TH1F* IsoRectangle_[ncell][nlayer];
  TH1F* IsoTriangleUpper_[ncell][nlayer];
  TH1F* IsoTriangleLower_[ncell][nlayer];
  
  TH2F* SeedVsIsolation_upper_[ncell][nlayer];
  TH2F* SeedVsIsolation_lower_[ncell][nlayer];
  TH2F* SeedVsIsolation_square_[ncell][nlayer];

  TH1F* ADC_Signal_[SAMPLES][nlayer];
  TH1F* ADC_All_Signal_[nlayer];
  TH1F* ADC_Signal[nlayer];
  TH1F* ADC_Signal_AllLayer[SAMPLES];
  TH2F* ADC_vs_TS_Signal_[nlayer];
  
  TH1F* ADC_Ctrl_[SAMPLES][nlayer];
  TH1F* ADC_All_Ctrl_[nlayer];
  TH1F* ADC_Ctrl[nlayer];
  TH1F* ADC_Ctrl_AllLayer[SAMPLES];
  TH2F* ADC_vs_TS_Ctrl_[nlayer];

  std::string  nameDetector_, DigiSource_;
  HGCalDDDConstants     *hgcons_;
  bool debug;
  explicit DigiAnalyzer(const edm::ParameterSet&);
  ~DigiAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  //  std::pair<HGCEEDetId, int> Max_Element(ADCInLayer adclayer);
  
  template<class T1, class T2> 
  void HGCDigiSaver(T1 detId, ForwardSubdetector subdet, const HGCalGeometry& geom0, const T2 it);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void DeclareHistograms();
  
  bool isSignal, isCtrl;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  // ----------member data ---------------------------
};
