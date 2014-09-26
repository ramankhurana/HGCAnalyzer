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

#include "DigiAnalyzer.h"

DigiAnalyzer::DigiAnalyzer(const edm::ParameterSet& iConfig)
{
  nameDetector_  = iConfig.getParameter<std::string>("DetectorName");
  DigiSource_ = iConfig.getParameter<std::string>("DigiSource");
  
  fs->make<TTree>("tree","tree");
  digiTree = fs->make<TTree>("digiTree","digiTree");
  DeclareHistograms();
  
  if (nameDetector_ == "HGCalEESensitive") {
    subdet = HGCEE;
  }
  else if (nameDetector_ == "HGCalHESiliconSensitive") subdet = HGCHEF;
  else subdet = HGCHEB;

  //now do what ever initialization is needed
  
  debug = false;
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
  
  isSignal = false;
  isCtrl   = false;
  for(int ilayer=0; ilayer<(int)AllADCInDet.size();ilayer++) AllADCInDet[ilayer].clear();
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag("genParticles"), genParticles);
  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[0] );
  
  edm::ESHandle<HGCalGeometry> geom;
  if (nameDetector_ == "HGCalEESensitive") {
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geom);
    
    if (!geom.isValid())
      std::cout << "Cannot get valid HGCalGeometry Object for " << "HGCalEESensitive" << std::endl;
    
    const HGCalGeometry& geom0 = *geom;
    edm::Handle<HGCEEDigiCollection> theHGCEEDigiContainers;
    iEvent.getByLabel("mix", DigiSource_, theHGCEEDigiContainers);
    
    for(HGCEEDigiCollection::const_iterator it =theHGCEEDigiContainers->begin();
	it !=theHGCEEDigiContainers->end(); ++it) {
      // skip hit if iterator  is empty 
      if(it->size()==0) continue;
      HGCEEDetId detId = it->id();
      
      if(it->sample(5).adc() ==0 ) continue;
      
      const GlobalPoint pos( std::move( geom0.getPosition( detId.rawId()) ) );
      
      float dR=deltaR(pos.eta(),pos.phi(),p.eta(),p.phi());
      isSignal = (dR<0.025);
      
      float idR=deltaR(pos.eta(),pos.phi(),-p.eta(),p.phi());
      isCtrl = (idR<0.025);
      
      //std::cout<<" dR = "<<dR<<std::endl;
      if(isSignal){
	int ilayer = detId.layer() - 1;
	std::pair<HGCEEDetId, int> tmpPair ( detId, it->sample(5).adc() );
	AllADCInDet[ilayer].insert(std::pair<HGCEEDetId, int>(tmpPair));
	std::cout<<" detId = "<<detId		 <<" adc = "<<it->sample(5).adc()		 <<std::endl;
      }
      
      
                  
      
      if(!isSignal && !isCtrl) continue;
      HGCDigiSaver(detId, subdet, geom0, it);
      
    } // end of digi collection
      
    for (int ilayer=0; ilayer<31; ilayer++) {
      maximumCell_ = 0;
      std::cout<<" size of each map for layer = "<<ilayer<<"   "<<AllADCInDet[ilayer].size()<<std::endl;
      ADCInLayer LayerADC = AllADCInDet[ilayer]; 
      std::pair<HGCEEDetId, int> maximumDetId (isocomputation->Max_Element(LayerADC) ) ;
      maximumCell_ = maximumDetId.second;
      SeedCell_[ilayer]->Fill(maximumCell_);
      std::cout<<" maximumCell_ in this layer "<<ilayer<<" is = "<<maximumCell_<<std::endl;
      std::vector<int> adciso = isocomputation->IsolatationRectangle(hgctopology_,LayerADC,ncell,ncell,1);
      for (int icell=0; icell<(int)adciso.size(); icell++) {
	IsoRectangle_[icell][ilayer]->Fill(adciso[icell]);
	SeedVsIsolation_square_[icell][ilayer]->Fill(maximumCell_,adciso[icell]);
      }
      
      std::vector<std::pair<int,int>> triangle_iso(isocomputation->IsolatationTriangle(hgctopology_,LayerADC,ncell,ncell,1));
      for (int icell=0; icell<(int)triangle_iso.size(); icell++) {
	IsoTriangleUpper_[icell][ilayer]->Fill(triangle_iso[icell].first);
	IsoTriangleLower_[icell][ilayer]->Fill(triangle_iso[icell].second);
	SeedVsIsolation_upper_[icell][ilayer]->Fill(maximumCell_,triangle_iso[icell].first);
	SeedVsIsolation_lower_[icell][ilayer]->Fill(maximumCell_,triangle_iso[icell].second);
      }
      for (int i=0; i<(int)triangle_iso.size(); i++) std::cout <<"adciso = "
							       <<adciso[i]
							       <<" tringle_iso = "
							       <<"    "<<triangle_iso[i].first
							       <<"    "<<triangle_iso[i].second
							       <<"    "<<maximumCell_
							       <<std::endl;
      
    }
    
  }
  
  
  // -------------------------------------//
  //  --------- End of HGC EE ------------//
  // -------------------------------------//

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
      std::cout<<" sending infor of digi "<<std::endl;
      
      
      
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
  
  int layer=detId.layer();
  int ilayer = layer -1;
  
  for (int i=0; i<nSample; ++i) {
    HGCSample hgcSample = it->sample(i);
    uint16_t gain = hgcSample.gain();
    uint16_t adc = hgcSample.adc();
    double charge = adc*gain;
    
    if(isSignal) {
      ADC_Signal_AllLayer[i]->Fill((int)adc);
      ADC_Signal_[i][ilayer]->Fill((int)adc);
      ADC_All_Signal_[ilayer]->Fill(i,(int)adc);
      ADC_vs_TS_Signal_[ilayer]->Fill((int)adc,i);
      ADC_Signal[ilayer]->Fill((int)adc);
    }
    if(isCtrl){
      ADC_Ctrl_AllLayer[i]->Fill((int)adc);
      ADC_Ctrl_[i][ilayer]->Fill((int)adc);
      ADC_All_Ctrl_[ilayer]->Fill(i,(int)adc);
      ADC_vs_TS_Ctrl_[ilayer]->Fill((int)adc,i);
      ADC_Ctrl[ilayer]->Fill((int)adc);
    }
    digiCharge += charge;
  }

  
}



// ------------ method called once each job just before starting event loop  ------------
void 
DigiAnalyzer::beginJob()
{
  //digiTree->Branch("ADC_All", "std::vector<int>", &ADC_All);
  //digiTree->Branch("Sample", "std::vector<int>", &Sample);
  //digiTree->Branch("ADC_TS1", "std::vector<int>", &ADC_TS1);
  //digiTree->Branch("ADC_TS2", "std::vector<int>", &ADC_TS2);
  //digiTree->Branch("ADC_TS3", "std::vector<int>", &ADC_TS3);
  //digiTree->Branch("Gain","std::vector<int>",&Gain);
  
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

void 
DigiAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESTransientHandle<DDCompactView> pDD;
  iSetup.get<IdealGeometryRecord>().get( pDD );
  const DDCompactView & cview = *pDD;
  hgcons_ = new HGCalDDDConstants(cview, nameDetector_);
  hgctopology_ = new HGCalTopology(*hgcons_,subdet,0);

}


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


void DigiAnalyzer::DeclareHistograms(){
  
  for(int i=0; i<ncell; i++){
    for(int ilayer=0; ilayer<nlayer; ilayer++){
    TString postname = Form("_ncell_%d_Layer_%d",i,ilayer);
    TString title    = Form(" %d cell %d;ADC;# of events (a.u.) ",ilayer+1,i+1);
    IsoRectangle_[i][ilayer] = fs->make<TH1F> ("Iso_Rectangle_"+postname, "Iso^{rectangle} : Layer"+title, 50,0,50);
    IsoTriangleUpper_[i][ilayer] = fs->make<TH1F> ("IsoTriangleUpper_"+postname, "Iso^{triangle}_{upper} : Layer"+title, 50,0,50);
    IsoTriangleLower_[i][ilayer] = fs->make<TH1F> ("IsoTriangleLower_"+postname, "Iso^{triangle}_{lower} : Layer"+title, 50,0,50);
    SeedVsIsolation_upper_[i][ilayer]  = fs->make<TH2F> ("SeedVsIsolation_upper_"+postname,"SeedVsIsolation_;Seed;Isolation",50,0,50, 50,0,50 );
    SeedVsIsolation_lower_[i][ilayer]  = fs->make<TH2F> ("SeedVsIsolation_lower_"+postname,"SeedVsIsolation_;Seed;Isolation",50,0,50, 50,0,50 );
    SeedVsIsolation_square_[i][ilayer] = fs->make<TH2F> ("SeedVsIsolation_square_"+postname,"SeedVsIsolation_;Seed;Isolation",50,0,50, 50,0,50 );
    }
  }

  for(int i=0; i<SAMPLES; i++){
    TString postname = Form("_TS_%d_",i);
    TString title    = Form(" %d;ADC;# of events (a.u.) ",i+1);
    ADC_Signal_AllLayer[i] = fs->make<TH1F> ("ADC_Signal_AllLayer"+postname,"ADC : TS"+title,50,0,50);
    ADC_Ctrl_AllLayer[i] = fs->make<TH1F> ("ADC_Ctrl_AllLayer"+postname,"ADC : TS"+title,50,0,50);
  }
  for(int ilayer=0; ilayer<nlayer; ilayer++){
    for(int i=0; i<SAMPLES; i++){
      TString postname = Form("_TS_%d_Layer_%d",i,ilayer);
      TString title    = Form(" %d TS %d;ADC;# of events (a.u.) ",ilayer+1,i+1);
      ADC_Signal_[i][ilayer] = fs->make<TH1F> ("ADC_Signal_"+postname, "ADC : Layer"+title, 50,0,50);
      ADC_Ctrl_[i][ilayer]   = fs->make<TH1F> ("ADC_Ctrl_"+postname, "ADC : Layer"+title,50,0,50);
    }
    TString postname1 = Form("_Layer_%d",ilayer);
    TString title1    = Form(" %d;ADC;# of events (a.u.) ",ilayer+1);
    SeedCell_[ilayer]    = fs->make<TH1F> ("SeedCell_"+postname1,"ADC counts in seed cell : Layer"+title1,50,0,50);
    ADC_Signal[ilayer]   = fs->make<TH1F> ("ADC_Signal_"+postname1,"ADC : Layer"+title1,50,0,50);
    ADC_Ctrl[ilayer]     = fs->make<TH1F> ("ADC_Ctrl_"+postname1,"ADC : Layer"+title1, 50,0,50);
    
    title1 = Form(" %d;ADC; TS",ilayer+1);
    ADC_vs_TS_Signal_[ilayer] = fs->make<TH2F> ("ADC_vs_TS_Signal_"+postname1,"ADC : Layer"+title1,100,-0.5,99.5, 6, 0., 6.);
    ADC_vs_TS_Ctrl_[ilayer]   = fs->make<TH2F> ("ADC_vs_TS_Ctrl_"+postname1,"ADC : Layer"+title1,100,-0.5,99.5, 6, 0., 6.);
    
    title1 = Form(" %d;TS; ADC ",ilayer+1);
    ADC_All_Signal_[ilayer]   = fs->make<TH1F> ("ADC_All_Signal_"+postname1,"ADC : Layer"+title1, 6,0,6);
    ADC_All_Ctrl_[ilayer]     = fs->make<TH1F> ("ADC_All_Ctrl_"+postname1,"ADC : Layer"+title1, 6,0,6);
  }
}

/*
std::pair<HGCEEDetId, int> DigiAnalyzer::Max_Element (std::map<HGCEEDetId, int> map_){
  
  std::cout<<" size of map = "<<map_.size()<<std::endl;
  std::map<HGCEEDetId,int>::const_iterator iter = map_.begin();
  int maximum_= iter->second;
  HGCEEDetId   key_;
  for(iter=map_.begin() ; iter != map_.end(); ++iter){
    //std::cout<< " finding maximum = "<< key_<<"  "<<maximum_<<std::endl;
    if(iter->second > maximum_){
      maximum_ = iter->second;
      key_     = iter->first;
    }
  }
  return std::pair<HGCEEDetId, int>(key_,maximum_);
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(DigiAnalyzer);


