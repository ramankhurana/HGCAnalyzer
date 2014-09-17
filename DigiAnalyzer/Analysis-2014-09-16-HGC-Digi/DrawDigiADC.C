void DrawDigiADC(){
  bool signal = true;
  bool PU     = false;
  /*
  // Trial using Keys :: Still under development
  
  TFile* f = new TFile("DigiTree.root","READ");
  
  TIter nextkey(f->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
  TString dirname = key->GetName();
  std::cout<<" isfolder = "<<key->IsFolder()<<std::endl;
  if(key->IsFolder()) {
  std::cout<<" dir =  "<<key->GetName()<<std::endl;  
  TDirectoryFile* fdir = (TDirectoryFile*) f->Get(key->GetName());
  TIter nextkeyOfdir(fdir->GetListOfKeys());
  TKey* keydir;
  while (keydir = (TKey*)nextkeyOfdir){
  std::cout<<" key name = "<<std::endl;
  
  }//while (keydir = (TKey*)nextkeyOfdir){
  }
  }
  
  */
  
  //TString files[]={"DigiTree.root","DigiTree.root","DigiTree.root","DigiTree.root"};
  
  //TString files[]={"DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0_140PU__TauValue-0.0.root.root", "DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0_140PU__TauValue-5.0.root.root", 
  //		  "DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0_140PU__TauValue-10.0.root.root", "DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0_140PU__TauValue-20.0.root.root"};

//  TString files[]={"DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0_140PU__TauValue-0.0.root.root", 
//		   "DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0_140PU__TauValue-5.0.root.root", 
//		   "DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0_140PU__TauValue-10.0.root.root",
//		   "DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0_140PU__TauValue-20.0.root.root"};
//  
  //TString files[]={"DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0__TauValue-0.0.root.root", "DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0__TauValue-5.0.root.root", 
  //		   "DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0__TauValue-10.0.root.root", "DigiHisto-gen-sim-digi-Pt-10.0_Eta-2.0__TauValue-20.0.root.root"};

  TString files[]={"DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0__TauValue-0.0.root.root", 
		   "DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0__TauValue-5.0.root.root", 
		   "DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0__TauValue-10.0.root.root", 
		   "DigiHisto-gen-sim-digi-SingleEle-Pt-10.0_Eta-2.0__TauValue-20.0.root.root"};
  

  TString legend[]={"#tau=0","#tau=5","#tau=10","#tau=20"};
  Int_t nevents;
  if(PU) nevents = 150;
  else nevents = 300;
  
  Int_t nlayers=30;
  Int_t colors[]={1,2,4,5};
  Int_t nfiles=sizeof(files) / sizeof(files[0]) ;
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  
  std::vector<TString> histname;
  histname.clear();
  for(int i=0; i<6; i++){
    if(signal) histname.push_back(Form("ADC_Signal_AllLayer_TS_%d_",i));
    else histname.push_back(Form("ADC_Ctrl_AllLayer_TS_%d_",i));
  }
  
  for(int i=0; i<30; i++) {
    if(signal){
      histname.push_back(Form("ADC_All_Signal__Layer_%d",i));
      histname.push_back(Form("ADC_Signal__Layer_%d",i));
    }
    else{
      histname.push_back(Form("ADC_All_Ctrl__Layer_%d",i));
      histname.push_back(Form("ADC_Ctrl__Layer_%d",i));
    }
  }
  
  for(int i=0; i<30; i++) {
    for(int j=0; j<6; j++){
      if(signal) histname.push_back(Form("ADC_Signal__TS_%d_Layer_%d",j,i));
      else histname.push_back(Form("ADC_Ctrl__TS_%d_Layer_%d",j,i));
      
    }
  }

  std::cout<<" histname size = "<<histname.size()<<std::endl;;
  
  TFile* fileIn;
  for(int i=0; i<histname.size();i++){
      
    TCanvas* c = new TCanvas("canvas", "canvas",286,86,600,600);
    //TCanvas* c = new TCanvas("canvas", "canvas",0,0,1300,700);
    c->Range(-68.75,-7.5,856.25,42.5);
    c->SetLogy(1);
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
    

    TLegend* leg;
    leg = new TLegend(0.78,0.7,0.97,0.9,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(22);
    
    TString pdfname;
    for(int ifile=0; ifile<nfiles;ifile++ ){
      fileIn = new TFile(files[ifile],"open");
      TH1F *h1 = (TH1F*) fileIn->Get("digianalyzerEE/"+histname[i]);
      
      h1->GetXaxis()->SetTitleFont(22);
      h1->GetYaxis()->SetTitleFont(22);
      h1->GetXaxis()->SetLabelFont(22);
      h1->GetYaxis()->SetLabelFont(22);
      h1->GetXaxis()->SetTitleSize(0.045);
      h1->GetYaxis()->SetTitleSize(0.045);
      h1->GetXaxis()->SetLabelSize(0.045);
      h1->GetYaxis()->SetLabelSize(0.045);
      h1->SetLineColor(colors[ifile]);
      h1->SetLineWidth(3);
      
      h1->SetTitleFont(22);
      double max = h1->GetMaximum();
      //h1->SetMaximum(1.2*max);
      
      leg->AddEntry(h1,legend[ifile],"F");
      pdfname = h1->GetName();
      if(i<6) h1->Scale(1.0/(nevents*nlayers));
      if(ifile==0) h1->Draw();
      else h1->Draw("same");
    }
    leg->Draw();
    if(PU){
      c->SaveAs("DigiPlots140PUEle/"+pdfname+".png");
      c->SaveAs("DigiPlots140PUEle/"+pdfname+".pdf");
    }
    else {
      c->SaveAs("DigiPlotsEle/"+pdfname+".png");
      c->SaveAs("DigiPlotsEle/"+pdfname+".pdf");
    }

    if(fileIn->IsOpen()) fileIn->Close();
  }
}
