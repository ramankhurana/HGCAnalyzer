#define DetectorEffects_cxx
#include "../interface/DetectorEffects.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
void DetectorEffects::Loop(TString outfilename)
{
  using namespace std;
  const double m_pi = 3.1415926535 ;
  OutputFileName = outfilename;
  f = new TFile(OutputFileName,"RECREATE");
  outTree_ = new TTree("outTree_","outTree_");
  
  MakeBranches();
  MakeHistos();
  
  std::cout<<" begin loop"<<std::endl;
  bool debug = false;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Clear();
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(debug) std::cout<< "n muons = "<<genjetngenMuons<<std::endl;
    if( genjetngenMuons ==2 ) {

      //compute dPt before filling the histograms 
      Float_t tmpsumdPt=0;
      for (auto i = 0;i < AK5nJet; i++){
	if((*AK5genjetEn)[i]<-900) continue; // don't use jet when there is no matched genjet for a given PFJet
	RECOJET_p4.SetPtEtaPhiE( (*AK5jetPt)[i],
				 (*AK5jetEta)[i],
				 (*AK5jetPhi)[i],
				 (*AK5jetEn)[i] );
	GENJET_p4.SetPxPyPzE( (*AK5genjetPx)[i],
			      (*AK5genjetPy)[i],
			      (*AK5genjetPz)[i],
			      (*AK5genjetEn)[i]);
	
	Float_t tmpDeltaPt_ = (GENJET_p4.Pt()-RECOJET_p4.Pt());
	tmpsumdPt =+ (tmpDeltaPt_) ;
      }//compute dPt before filling the histograms ends here 

      std::cout<<" value of sum dpt = "<<tmpsumdPt<<std::endl;

      Float_t HT=0;
      Float_t sumdPt =0.;
      for (auto j = 0;j < AK5nJet; j++){
	if((*AK5genjetEn)[j]<-900) continue; // don't use jet when there is no matched genjet for a given PFJet
	TLorentzVector RECOJET_p4;
	RECOJET_p4.SetPtEtaPhiE( (*AK5jetPt)[j],
				 (*AK5jetEta)[j],
				 (*AK5jetPhi)[j],
				 (*AK5jetEn)[j] );
	
	
	
	TLorentzVector GENJET_p4;
	GENJET_p4.SetPxPyPzE( (*AK5genjetPx)[j],
			      (*AK5genjetPy)[j],
			      (*AK5genjetPz)[j],
			      (*AK5genjetEn)[j]);
	
	DeltaPt_ = (GENJET_p4.Pt()-RECOJET_p4.Pt())/GENJET_p4.Pt();
	
	HT =+ RECOJET_p4.Pt(); // compute sum pT of all the jets
	sumdPt =+ (DeltaPt_*GENJET_p4.Pt()) ;
	
	//Fill branches to New Tree;
	dpT_.push_back((double)DeltaPt_*GENJET_p4.Pt()) ;
	dpT_Over_pT_.push_back((double)DeltaPt_);
	GenJetpT_.push_back((double)GENJET_p4.Pt());
	GenJeteta_.push_back((double)GENJET_p4.Eta());
	GenJetphi_.push_back((double)GENJET_p4.Phi());
	
	
	// This will make sure that only outliers will be filled in the histograms
	//if( tmpsumdPt > 250.0 && pfMetRawPt > 250.0) {
	//if( (DeltaPt_*GENJET_p4.Pt()) > 200.0 && pfMetRawPt > 200.0) {
	delta->Fill(DeltaPt_);
	dpT->Fill((DeltaPt_*GENJET_p4.Pt()));
	
	deltapt_vs_eta->Fill( GENJET_p4.Eta(),DeltaPt_);
	deltapt_vs_phi->Fill( GENJET_p4.Phi(),DeltaPt_);
	//eta_vs_phi_profile_dPt->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), (GENJET_p4.Energy() - RECOJET_p4.Energy()));
	eta_vs_phi_profile_dPt->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), DeltaPt_);
	MET_vs_dPt->Fill(pfMetRawPt,(DeltaPt_*GENJET_p4.Pt()));
	dpT_vs_genJetpT->Fill(DeltaPt_*GENJET_p4.Pt(),GENJET_p4.Pt());
	dpT_vs_recoJetpT->Fill(DeltaPt_*GENJET_p4.Pt(),RECOJET_p4.Pt());
	
	GenEM_vs_RecoEM->Fill((*AK5genjetEM)[j]/GENJET_p4.Energy(), (*AK5jetPhoEF)[j]);
	GenHAD_vs_RecoChHAD->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j]);
	GenHAD_vs_RecoHAD->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j] + (*AK5jetNHadEF)[j]);
	

	// Dilute with dpT > 200 GeV && MET > 200 GeV
	if( (DeltaPt_*GENJET_p4.Pt()) > 200.0 && pfMetRawPt > 200.0) { 
	  eta_vs_phi_profile_dPt_DilutedpT->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), DeltaPt_);
       	}
	
	//Dilute with dphi (MET, jet pT) < 0.5
	double phi1 = GENJET_p4.Phi();
	double phi2 = pfMetRawPhi;;
	double result = phi1 - phi2;
	if(result>m_pi) result -= 2*m_pi;
	else if (result <= -m_pi) result += 2*m_pi ;
	else result = result;
	dPhi_MET_Jet->Fill(result);
	dPhi_vs_MET->Fill(result, pfMetRawPt);
	if(result <0.5){
	  eta_vs_phi_profile_dPt_Dilute_dphi->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), DeltaPt_);
	  MET_vs_dPt_Dilute_dphi->Fill(pfMetRawPt,(DeltaPt_*GENJET_p4.Pt()));
	}
	
      }//for (auto j = 0;j < AK5nJet; j++){
      //if( tmpsumdPt > 350.0 && pfMetRawPt > 350.0) {
      //if( (DeltaPt_*GENJET_p4.Pt()) > 200.0 && pfMetRawPt > 200.0) {
      MET_vs_HT->Fill(pfMetRawPt,HT);
      MET_vs_SumdPt->Fill(pfMetRawPt,sumdPt);
      SumdpT->Fill(sumdPt);
      //      }
      
      // Fill New Branches;
      MET_ = (float) pfMetRawPt;
      HT_ = (float) HT;
      SumdpT_ = (float) sumdPt;
      
      Run_ = info_runId;
      Lumi_ = info_lumiSection;
      Event_ = info_eventId;
      outTree_->Fill(); 
    }    // exact two gen muons
    
    //electrons.pVector = (TLorentzVector*) RK_Electron_4Momentum->At(i);
    if(jentry%100000==1) std::cout<<" event number  = "<<jentry<<std::endl;
  }
  
  
  f->cd();
  //outTree_->Write();
  
  delta->Write();
 
  SumdpT->Write();
  
  dpT->Write();
  deltapt_vs_eta->Write();
  deltapt_vs_phi->Write();
  eta_vs_phi_profile_dPt->Write();
  MET_vs_HT->Write();
  MET_vs_SumdPt->Write();
  MET_vs_dPt->Write();
  eta_vs_phi_profile_dPt_DilutedpT->Write();
  eta_vs_phi_profile_dPt_Dilute_dphi->Write();
  MET_vs_dPt_Dilute_dphi->Write();
  dpT_vs_genJetpT->Write();
  dpT_vs_recoJetpT->Write();
  GenEM_vs_RecoEM->Write();
  GenHAD_vs_RecoChHAD->Write();
  GenHAD_vs_RecoHAD->Write();
  dPhi_MET_Jet->Write();
  dPhi_vs_MET->Write();
  f->Close();
  std::cout<<" file closed and job finished"<<std::endl;
}




void DetectorEffects::MakeBranches(){
  
  outTree_->Branch("MET_",&MET_,"MET_/F");
  outTree_->Branch("HT_",&HT_,"HT_/F");
  outTree_->Branch("SumdpT_",&SumdpT_,"SumdpT_/F");
  
  outTree_->Branch("dpT_","std::vector<double>",&dpT_);
  outTree_->Branch("dpT_Over_pT_","std::vector<double>",&dpT_Over_pT_);
  outTree_->Branch("GenJetpT_","std::vector<double>",&GenJetpT_);
  outTree_->Branch("GenJeteta_","std::vector<double>",&GenJeteta_);
  outTree_->Branch("GenJetphi_","std::vector<double>",&GenJetphi_);
  outTree_->Branch("Run_",&Run_,"Run_/I");
  outTree_->Branch("Lumi_",&Lumi_,"Lumi_/I");
  outTree_->Branch("Event_",&Event_,"Event_/I");
}


void DetectorEffects::MakeHistos(){
  delta = new TH1F("delta","delta;(p_{T}^{genJet}-p_{T}^{recoJet})/p_{T}^{genJet};# of entries",120,-3,3.);
  SumdpT = new TH1F("SumdpT","SumdpT;#Sigma #Delta p_{T};# of Events",1000,-1000,1000);
  dpT = new TH1F("dpT","dpT; #Delta p_{T};# of Events",1000,-1000,1000);
  dPhi_MET_Jet = new TH1F("dPhi_MET_Jet","dPhi_MET_Jet;#Delta #phi (MET,jet);# of Events", 100,-20,20);
  
  deltapt_vs_eta = new TH2F("deltapt_vs_eta","deltapt_vs_eta; #eta_{genJet} ; #Delta p_{T}/p_{T}^{genJet}",100,-3.14,3.14, 400,-10.,10.);
  deltapt_vs_phi = new TH2F("deltapt_vs_phi","deltapt_vs_phi; #phi_{genJet}; #Delta p_{T}/p_{T}^{genJet} ",100,-3.14,3.14, 400,-10.,10.);
  eta_vs_phi_profile_dPt = new TProfile2D("eta_vs_phi_profile_dPt","eta_vs_phi_profile_dPt;#eta^{genJet};#phi_{genJet}",100,-3.14,3.14, 100,-3.14,3.14); 
  MET_vs_HT = new TH2F("MET_vs_HT","MET_vs_HT;MET; H_{T}",200,0,1000,200,0,1000);
  MET_vs_SumdPt = new TH2F("MET_vs_SumdPt","MET_vs_SumdPt;MET; #Sigma #Delta p_{T}",200,0,1000,200,0,1000);  
  MET_vs_dPt    = new TH2F("MET_vs_dPt","MET_vs_dPt;MET;#Delta p_{T}",200,0,1000,200,0,1000);
  
  dpT_vs_genJetpT = new TH2F("dpT_vs_genJetpT","dpT_vs_genJetpT",1000,-1000,1000, 1000,-1000,1000);
  dpT_vs_recoJetpT = new TH2F("dpT_vs_recoJetpT","dpT_vs_recoJetpT",1000,-1000,1000, 1000,-1000,1000);
  
  eta_vs_phi_profile_dPt_DilutedpT = new TProfile2D("eta_vs_phi_profile_dPt_DilutedpT","eta_vs_phi_profile_dPt_DilutedpT;#eta^{genJet};#phi_{genJet}",100,-3.14,3.14, 100,-3.14,3.14); 
  eta_vs_phi_profile_dPt_Dilute_dphi = new TProfile2D("eta_vs_phi_profile_dPt_Dilute_dphi","eta_vs_phi_profile_dPt_Dilute_dphi;#eta^{genJet};#phi_{genJet}",100,-3.14,3.14, 100,-3.14,3.14); 
  MET_vs_dPt_Dilute_dphi    = new TH2F("MET_vs_dPt_Dilute_dphi","MET_vs_dPt_Dilute_dphi;MET;#Delta p_{T}",200,0,1000,200,0,1000);
  
  GenEM_vs_RecoEM = new TH2F("GenEM_vs_RecoEM","GenEM_vs_RecoEM;EM_{gen};EM_{reco}",100,0,2,100,0,2);
  GenHAD_vs_RecoChHAD = new TH2F("GenHAD_vs_RecoChHAD","GenHAD_vs_RecoChHAD;Had_{gen};ChHad_{reco}",100,0,2,100,0,2);
  GenHAD_vs_RecoHAD = new TH2F("GenHAD_vs_RecoHAD","GenHAD_vs_RecoHAD;Had_{gen};ChHad_{reco}",100,0,2,100,0,2);
  
  dPhi_vs_MET = new TH2F("dPhi_vs_MET","dPhi_vs_MET;#Delta #phi; MET",100,-20,20,500,0,1000);
}



void DetectorEffects::Clear(){
  MET_ = 0.;
  HT_ = 0.;
  SumdpT_ = 0.;
  dpT_.clear();
  dpT_Over_pT_.clear();
  GenJetpT_.clear();
  GenJeteta_.clear();
  GenJetphi_.clear();
  Run_ = 0;
  Lumi_ = 0;
  Event_ = 0;
}
