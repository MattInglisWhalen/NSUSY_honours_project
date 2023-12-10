#ifndef LOOK_AT_HEP_EVENTS_H
#define LOOK_AT_HEP_EVENTS_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::look_at_HEP_events(){

  TFile* f = new TFile("HEP_events.root","RECREATE");

  // pre-cut histograms

  TH1I* h_nphot = new TH1I("h_nphot","Number of Photons",500,0,500);
  TH1F* h_photpt = new TH1F("h_photpt","p_{T} Distributions of Photons",50,0,100);

  TH1I* h_nelec = new TH1I("h_nelec","Number of Electrons",10,0,10);
  TH1F* h_elecpt = new TH1F("h_elecpt","p_{T} Distribution of Electrons",100,0,200);

  TH1I* h_nmuon = new TH1I("h_nmuon","Number of Muons",10,0,10);
  TH1F* h_muonpt = new TH1F("h_muonpt","p_{T} Distributions of Muons",100,0,300);

  TH1I* h_ntau = new TH1I("h_ntau","Number of Tauons",10,0,10);
  TH1F* h_taupt = new TH1F("h_taupt","p_{T} Distributions of Tauons",100,0,100);

  TH1I* h_nneut = new TH1I("h_nneut","Number of Neutrinos",10,0,10);
  TH1F* h_neutpt = new TH1F("h_neutpt","p_{T} Distributions of Neutrinos",100,0,300);

  TH1I* h_njet_uncut  = new TH1I("h_njet_uncut","Number of Jets before Cuts",50,0,50);
  TH1F* h_jetpt_uncut = new TH1F("h_jetpt_uncut","Jet p_{T} before Cuts",100,0,500);
  TH1F* h_jeteta_uncut = new TH1F("h_jeteta_uncut","Jet #eta before Cuts",100,-12.0,12.0);
  TH1F* h_massj_uncut = new TH1F("h_massj_uncut","Jet's Mass Before Cuts",100,0,200);

  TH1F* h_met_uncut = new TH1F("h_metpt_uncut","E_{T}^{miss} Before Cuts",100,0,300);

  // post-cut histograms

  TH1I* h_njet_cut  = new TH1I("h_njet_cut","Number of Jets after Cuts",50,0,50);
  TH1F* h_jetpt_cut = new TH1F("h_jetpt_cut","Jet p_{T} after Cuts",100,0,500);
  TH1F* h_jeteta_cut = new TH1F("h_jeteta_cut","Jet #eta after Cuts",100,-12.0,12.0);
  TH1F* h_massj_cut = new TH1F("h_massj_cut","Jet's Mass After Cuts",100,0,200);

  TH1F* h_met_cut = new TH1F("h_metpt_cut","E_{T}^{miss} After Cuts",100,0,300);

  TH1F* h_1jincl80 = new TH1F("h_1jincl80","#sigma_{1} vs. p_{T}^{cut}",80,0,80);
  TH1F* h_1jincl500 = new TH1F("h_1jincl500","#sigma_{1} vs. p_{T}^{cut}",250,0,500);

  TH1F* h_ptveto80 = new TH1F("h_ptveto80","#sigma_{0} vs. p_{T}^{cut}",80,0,80);
  TH1F* h_ptveto500 = new TH1F("h_ptveto500","#sigma_{0} vs. p_{T}^{cut}",250,0,500);

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    progress_bar(ientry,nEvents);
    read_event(ientry);

    // Fill histograms with uncut data

    h_nphot->Fill(nphot);
    for(Int_t i=0;i<nphot;++i){
      h_photpt->Fill( sphot[i].pt() );
    }

    h_nelec->Fill(nelec);
    for(Int_t i=0;i<nelec;++i){
      h_elecpt->Fill( selec[i].pt() );
    }

    h_nmuon->Fill(nmuon);
    for(Int_t i=0;i<nmuon;++i){
      h_muonpt->Fill( smuon[i].pt() );
    }

    h_ntau->Fill(ntau);
    for(Int_t i=0;i<ntau;++i){
      h_taupt->Fill( stau[i].pt() );
    }

    h_nneut->Fill(nneut);
    for(Int_t i=0;i<nneut;++i){
      h_neutpt->Fill( sneut[i].pt() );
    }
    
    h_njet_uncut->Fill(njet);
    for(Int_t i=0;i<njet;++i){
      h_jetpt_uncut->Fill( sjet[i].pt() );
      h_jeteta_uncut->Fill( sjet[i].eta() );
      h_massj_uncut->Fill( sjet[i].mass() );
    }
   
    h_met_uncut->Fill( smissing.pt() );

    Int_t vetoflag;
    for(Int_t i=0;i<80;i+=1){
      vetoflag=0;
      if(sjet[0].pt()>i){
        h_1jincl80->Fill(i);
      }
      for(Int_t j=0;j<njet;++j){
        if( (sjet[j].pt()>i) && (sjet[j].abseta()<4.5) ){
          vetoflag=1;
        }
      }
      if(!vetoflag){h_ptveto80->Fill(i);}
    }
    for(Int_t i=1;i<500;i+=2){
      vetoflag=0;
      if(sjet[0].pt()>i){
        h_1jincl500->Fill(i);
      }
      for(Int_t j=0;j<njet;++j){
        if( (sjet[j].pt()>i) && (sjet[j].abseta()<4.5) ){
          vetoflag=1;
        }
      }
      if(!vetoflag){h_ptveto500->Fill(i);}
    }

    // Perform the kinematic cuts on the final state particles
    // Require jets to have pt>20 and |eta|<2.5

    //reconstruct(sphot,nphot,5.,2.5);
    //reconstruct(selec,nelec,10.,2.5);
    //reconstruct(smuon,nmuon,10.,2.5);
    //reconstruct(stau,ntau,10.,2.5);
    //reconstruct(sneut,nneut,5.,2.5);
    reconstruct(sjet,njet,20.,2.5);
      
    // Fill histograms with cut data

    h_njet_cut->Fill(njet);
    for(Int_t i=0;i<njet;++i){
      h_jetpt_cut->Fill( sjet[i].pt() );
      h_jeteta_cut->Fill( sjet[i].eta() );
      h_massj_cut->Fill( sjet[i].mass() );
    }
    recalculate_smissing();

    h_met_cut->Fill( smissing.pt() );

  }

  h_1jincl80->Scale(24.99/h_1jincl80->GetBinContent(1));
  h_1jincl500->Scale(24.99/h_1jincl500->GetBinContent(1));

  h_ptveto80->Scale(24.99/h_ptveto500->GetBinContent(250));
  h_ptveto500->Scale(24.99/h_ptveto500->GetBinContent(250));

  f->Write();

};

#endif
