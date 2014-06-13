// -*- C++ -*-
//
// Package:    HGCSimHitsAnalyzer
// Class:      HGCSimHitsAnalyzer
// 
/**\class HGCSimHitsAnalyzer HGCSimHitsAnalyzer.cc Validaion/HGCalValidation/plugins/HGCSimHitsAnalyzer.cc

 Description: Validates SimHits of High Granularity Calorimeter

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raman Khurana
////      and Kalyanmoy Chatterjee
//         Created:  Fri, 31 Jan 2014 18:35:18 GMT
// $Id$
//
//

// system include files
#include "HGCSimHitsAnalyzer.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DetectorDescription/Core/interface/DDExpandedView.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "TVector3.h"
#include <cmath>
#include "CommonTools/UtilAlgos/interface/TFileService.h"



HGCSimHitsAnalyzer::HGCSimHitsAnalyzer(const edm::ParameterSet& iConfig){
  //now do what ever initialization is needed
  //dbe_           = edm::Service<DQMStore>().operator->();
  nameDetector_   = iConfig.getParameter<std::string>("DetectorName");
  caloHitSource_  = iConfig.getParameter<std::string>("CaloHitSource");
  outputfilename_ = iConfig.getParameter<std::string>("OutFile");
  
  sampFracEE = 1.54338469862809234e+02 ;
  f = new TFile(outputfilename_.c_str(),"RECREATE");

  const char *name =  nameDetector_.c_str();
  tree = new TTree("hgcSimHitTree",name);

  //SimHitVector.reserve(2000);
}


HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer() {
 
}


// ------------ method called for each event  ------------
void HGCSimHitsAnalyzer::analyze(const edm::Event& iEvent, 
				 const edm::EventSetup& iSetup) {
  
  SimHitVector.clear();
  
  std::cout<<" inside analyze "<<std::endl;
  edm::Handle<edm::PCaloHitContainer> theCaloHitContainers;
  iEvent.getByLabel("g4SimHits", caloHitSource_, theCaloHitContainers);
  if (verbosity_>0) std::cout << " PcalohitItr = " 
			      << theCaloHitContainers->size() << std::endl;
  std::vector<PCaloHit>               caloHits;
  caloHits.insert(caloHits.end(), theCaloHitContainers->begin(), 
		  theCaloHitContainers->end());

  analyzeHits(caloHits);
  tree->Fill();
}

void HGCSimHitsAnalyzer::analyzeHits (std::vector<PCaloHit>& hits) {
  double energy_;       energy_     =  0.0;
  double mipInverse_;   mipInverse_ =  0.0;
  int nHit = hits.size();
  //for (int i=0; i<nHit; i++) {
  for (int i=0; i<nHit; i++) {
    if(i>1999) continue ;

    HGCalDetId detId   = HGCalDetId(hits[i].id());
    double energy      =   hits[i].energy();
    energy_ = energy_ + energy;
    
    double energyem    =   hits[i].energyEM();
    double energyhad   =   hits[i].energyHad();
    double time        =   hits[i].time();
    int    cell        =   detId.cell();
    int    subsector   =   detId.subsector();
    int    sector      =   detId.sector();
    int    layer       =   detId.layer();
    int    zside       =   detId.zside();
    
    //    if (verbosity_>1)
      std::cout << " energy = " << energy
				<< " energyem = "<<energyem
				<< " energyhad = " << energyhad
				<< " time = " << time
				<< " cell = " << cell
				<< " subsector = " << subsector
				<< " sector = " << sector
				<< " layer = " << layer
				<< " zside = " << zside << std::endl;
    std::vector<double> esumVector;
    energysum esum;
    if(time > 0){
      esum.e1000 += energy;
      esumVector.push_back(esum.e1000);
      if(time < 250) {
        esum.e250 += energy;
	esumVector.push_back(esum.e250);
        if(time < 100){
          esum.e100 += energy;
	  esumVector.push_back(esum.e100);
          if(time < 50){
            esum.e50 += energy;
	    esumVector.push_back(esum.e50);
            if(time < 25){
              esum.e25 += energy;
	      esumVector.push_back(esum.e25);
              if(time < 15){
                esum.e15 += energy;
		esumVector.push_back(esum.e15);
              }//if(time < 15)
            }//if(time < 25)
          }//if(time < 50)
        }//if(time < 100)
      } // if(time < 250)
    } // if(time < 1000)


    std::pair<float,float> xy = hgcons_->locateCell(cell,layer,subsector,false);
    const HepGeom::Point3D<float> lcoord(xy.first,xy.second,0);
    int      subs = (symmDet_ ? 0 : subsector);
    uint32_t id = HGCalDetId(ForwardEmpty,zside,layer,sector,subs,0).rawId();
    if (verbosity_>1) std::cout << HGCalDetId(id) << " in map " 
				<< (transMap_.find(id) != transMap_.end()) 
				<< std::endl;
    const HepGeom::Point3D<float> gcoord( transMap_[id]*lcoord );

    float globalx = gcoord.x();
    float globaly = gcoord.y();
    float globalz = gcoord.z();
    hitsinfo   hinfo;
    TVector3   HitGlobalCo;
    HitGlobalCo.SetXYZ(globalx, globalx, globalz);
    double velocity = 3*pow(10,11);
    double tof = HitGlobalCo.Mag()/velocity; 

    hinfo.x          = globalx;//xy.first;
    hinfo.y          = globaly;//xy.second;
    hinfo.z          = globalz;
    hinfo.time       = time - tof;
    hinfo.energy     = energy;
    hinfo.sector     = sector;
    hinfo.cell       = cell;
    hinfo.layer      = layer;
    hinfo.phi        = atan2(globaly,globalx);
    double theta     = acos(globalz/sqrt(globalx*globalx+globaly*globaly+globalz*globalz));
    hinfo.eta        = -log(tan(theta/2)); 
    hinfo.subdetname = nameDetector_;
    if (verbosity_>1) std::cout << " --------------------------   gx = " 
				<< globalx << " gy = "  << globaly   << " gz = "
				<< globalz << " phi = " << hinfo.phi << " eta = "
				<< hinfo.eta << std::endl;
    std::pair<hitsinfo,energysum> pair_tmp(hinfo,esum);
    SimHitVector.push_back(hinfo);

    std::cout<<" before fill histo "<<std::endl;
    for(unsigned int itimeslice = 0; itimeslice < esumVector.size(); itimeslice++ ) {
      FillHitsInfo(pair_tmp, itimeslice, esumVector.at(itimeslice));
    }
        
    mipInverse_ = mipInverse_ + OneOverMip(pair_tmp);
  }
  EnergyCalibrated_->Fill(energy_*sampFracEE);
}

void HGCSimHitsAnalyzer::FillHitsInfo(std::pair<hitsinfo,energysum> hits, unsigned int itimeslice, double esum){
  if (!geometrydefined_) return;
  
}


bool HGCSimHitsAnalyzer::defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH){

  const DDCompactView & cview = *ddViewH;
  hgcons_ = new HGCalDDDConstants(cview, nameDetector_);
  if (verbosity_>0) std::cout << "Initialize HGCalDDDConstants for " 
			      << nameDetector_ << " : " << hgcons_ <<std::endl;

  std::string attribute = "Volume"; 
  std::string value     = nameDetector_;
  DDValue val(attribute, value, 0);
  
  DDSpecificsFilter filter;
  filter.setCriteria(val, DDSpecificsFilter::equals);
  DDFilteredView fv(cview);
  fv.addFilter(filter);
  bool dodet = fv.firstChild();

  while (dodet) {
    const DDSolid & sol  = fv.logicalPart().solid();
    std::string name = sol.name();
    int isd = (name.find(nameDetector_) == std::string::npos) ? -1 : 1;
    if (isd > 0) {
      std::vector<int> copy = fv.copyNumbers();
      int nsiz = (int)(copy.size());
      int lay  = (nsiz > 0) ? copy[nsiz-1] : -1;
      int sec  = (nsiz > 1) ? copy[nsiz-2] : -1;
      int zp   = (nsiz > 3) ? copy[nsiz-4] : -1;
      if (zp !=1 ) zp = -1;
      const DDTrap & trp = static_cast<DDTrap>(sol);
      int subs = (trp.alpha1()>0 ? 1 : 0);
      symmDet_ = (trp.alpha1()==0 ? true : false);
      uint32_t id = HGCalDetId(ForwardEmpty,zp,lay,sec,subs,0).rawId();
      DD3Vector x, y, z;
      fv.rotation().GetComponents( x, y, z ) ;
      const CLHEP::HepRep3x3 rotation ( x.X(), y.X(), z.X(),
					x.Y(), y.Y(), z.Y(),
					x.Z(), y.Z(), z.Z() );
      const CLHEP::HepRotation hr ( rotation );
      const CLHEP::Hep3Vector h3v ( fv.translation().X(),
				    fv.translation().Y(),
				    fv.translation().Z()  ) ;
      const HepGeom::Transform3D ht3d (hr, h3v);
      transMap_.insert(std::make_pair(id,ht3d));
      if (verbosity_>1) std::cout << HGCalDetId(id) << " Transform using " 
				  << h3v << " and " << hr;
    }
    dodet = fv.next();
  }
  if (verbosity_>0) std::cout << "Finds " << transMap_.size() << " elements" 
			      << " and SymmDet_ = " << symmDet_ << std::endl;
  return true;
}


// ------------ method called once each job just before starting event loop  ------------
void HGCSimHitsAnalyzer::beginJob() {
  geometrydefined_ = false;
  symmDet_         = true;
  tree->Branch(("SimHitVector"),"std::vector<hitsinfo>",&SimHitVector);
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCSimHitsAnalyzer::endJob() {
  f->cd();
  EnergyCalibrated_->Write();
  tree->Write();
}

// ------------ method called when starting to processes a run  ------------

void HGCSimHitsAnalyzer::beginRun(edm::Run const& iRun, 
				     edm::EventSetup const& iSetup) {
  
  std::cout<<" inside beginrun"<<std::endl;
  if (!geometrydefined_) {
    edm::ESTransientHandle<DDCompactView> pDD;
    iSetup.get<IdealGeometryRecord>().get( pDD );
    geometrydefined_ = defineGeometry(pDD);
    
    std::cout<<" after geometrydefined_ = "<<geometrydefined_<<std::endl;
    unsigned int layers = hgcons_->layers(false);
    nlayers = layers;
    EnergyCalibrated_   = new TH1F("EnergyCalibrated_","Energy_",1000,0.0,1000.);
    
  }
  
}


// ------------ method called when ending the processing of a run  ------------
/*
void HGCSimHitsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void HGCSimHitsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void HGCSimHitsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCSimHitsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



float HGCSimHitsAnalyzer::OneOverMip(std::pair<hitsinfo,energysum> hit_){
  //float ETA = hit_.first.eta ;
  //float MIP      = 54.3/TMath::TanH(ETA); // in KeV
  float MIP     =  1.1;
  return( 1.0 / MIP); // 1000000 is because of GeV
}

#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
