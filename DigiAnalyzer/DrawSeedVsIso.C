#include "TH2F.h"
#include "TStyle.h"
#include <vector>
#include "TROOT.h"
#include <iostream>
#include <TCanvas.h>
#include "TFile.h"
using namespace std;
void DrawSeedVsIso(){
  TCanvas* c = new TCanvas("canvas", "canvas",286,86,600,600);
  c->Range(-68.75,-7.5,856.25,42.5);
  c->SetLogy(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.10);
  c->SetBottomMargin(0.20);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetPadTickX(3);
  gStyle->SetPadTickY(3);
  
  
  
  std::vector<TH2F*> HistoVector;
  HistoVector.clear();
  TFile* f = new TFile("DigiHisto-INPUTFILENAME.root","READ");
  std::vector<TString> midname;
  midname.clear();
  midname.push_back("square");
  midname.push_back("upper");
  midname.push_back("lower");

  TString prename ;
  TString postname;
  TString name    ;


  for(int ilayer = 0; ilayer < 30; ilayer++){
    for(int icell = 0; icell < 5; icell++){
      for(int imidname=0; imidname <(int)midname.size(); imidname++){
	prename  = "SeedVsIsolation_"+midname[imidname];
	postname = Form("__ncell_%d_Layer_%d",icell,ilayer);
	name     = prename + postname ;
	TH2F* h2 = (TH2F*) f->Get("digianalyzerEE/"+name);
	h2->Draw("colz");
	h2->SetTitleFont(22,"XYZ");
	h2->SetLabelFont(22,"XYZ");
	h2->SetTitleOffset(1.45,"Y");
	h2->SetTitleSize(0.045,"XYZ");
	h2->SetLabelSize(0.045,"XYZ");
	c->SaveAs("plots/"+name+".pdf");
	c->SaveAs("plots/"+name+".png");
	c->SaveAs("plots/"+name+".root");
      }
    }
  }
  f->Close();
}
