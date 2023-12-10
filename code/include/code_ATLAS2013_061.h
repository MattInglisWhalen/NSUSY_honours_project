#ifndef CODE_ATLAS2013_061_H
#define CODE_ATLAS2013_061_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_ATLAS2013_061(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnelec = new TH1I("hnelec","Number of Electrons",10,0,10);
  TH1I* hnmuon = new TH1I("hnmuon","Number of Muons",10,0,10);
  TH1I* hnlep  = new TH1I("hnlep","Number of Leptons",8,0,8);
  TH1I* hnrelep  = new TH1I("hnrelep","Number of Reconstructed Leptons",8,0,8);
  TH1I* hnjet  = new TH1I("hnjet","Number of Jets",25,0,25);
  TH1I* hnbtag = new TH1I("hnbtag","Number of total btags",15,0,15);
  TH1I* hntrack  = new TH1I("hntrack","Number of Tracks",300,0,300);
 
  TH1F* helecpt1 = new TH1F("helecpt1","P_t of First Electron",100,0,1000);
  TH1F* helecpt2 = new TH1F("helecpt2","P_t of Second Electron",100,0,1000);
  TH1F* hchargee = new TH1F("hchargee","Electron's Charge",20,-5,5);
  TH1F* hhademe = new TH1F("hhademe","Electron's Had/Em",100,0,0.25);

  TH1F* hmuonpt1 = new TH1F("hmuonpt1","P_t of First Muon",100,0,1000);
  TH1F* hmuonpt2 = new TH1F("hmuonpt2","P_t of Second Muon",100,0,1000);
  TH1F* hchargem = new TH1F("hhchargem","Muon's Charge",20,-5,5);

  TH1F* hjetpt1 = new TH1F("hjetpt1","P_t of First Jet",100,0,1000);
  TH1F* hjetpt2 = new TH1F("hjetpt2","P_t of Second Jet",100,0,1000);
  TH1F* hjetpt3 = new TH1F("hjetpt3","P_t of Third Jet",100,0,1000);
  TH1F* hjetpt4 = new TH1F("hjetpt4","P_t of Fourth Jet",100,0,1000);
  TH1F* hmassj = new TH1F("hmassj","Jet's Mass",100,0,200);
  TH1F* hhademj = new TH1F("hhademj","Jet's Had/Em",100,0,20);

  TH1F* hmetpt = new TH1F("hmetpt","Missing Transverse Energy",100,0,1000);
  TH1F* hmetphi= new TH1F("hmetphi","Phi of Missing Transverse Energy",32,-TMath::Pi(),TMath::Pi());  

  TH1F* hrmetpt = new TH1F("hrmetpt","Recalculated Missing Transverse Energy",100,0,1000);
  TH1F* hrmetphi= new TH1F("hrmetphi","Recalculated Phi of Missing Transverse Energy",32,-TMath::Pi(),TMath::Pi());  

  TH1F* hHT = new TH1F("hHT","HT",100,0,2500);
  TH1F* hHT_4j = new TH1F("hHT_4j","HT_4j",100,0,2500);
  TH1F* hmeff = new TH1F("hmeff","Effective Mass",100,0,3500);
  TH1F* hmeff_4j = new TH1F("hmeff_4j","Effective Mass with Leading 4 Jets",100,0,3500);
  TH1F* hdelphimin = new TH1F("hdelphimin","Minimum Phi Between a Jet and Missing Energy",100,0,TMath::Pi());
  TH1F* hET_meff_4j = new TH1F("hET_meff_4j","Missing Energy over 4-jet Meff",100,0,2);
  TH1F* hET_HT_4j = new TH1F("hET_HT_4j","Missing Energy over HT",100,0,50);

  TH1F* hHT_inc = new TH1F("hHT_inc","HT with lepton",100,0,2500);
  TH1F* hMT = new TH1F("hMT","Transverse Mass",100,0,2500);
  TH1F* hmeff_inc = new TH1F("hmeff_inc","Effective Mass with lepton",100,0,3500);
  TH1F* hET_HT_inc = new TH1F("hET_HT_inc","Missing Energy over HT with lepton",100,0,50);

  TH1I* hnjet50  = new TH1I("hnjet50","Number of Jets with pt>50",25,0,25);
  TH1I* hnbtag50 = new TH1I("hnbtag50","Number of total btags with pt>50",15,0,15);  

  TH1I* hcut=new TH1I("hcut","Where the cuts occur",20,0,20);
  TH1I* hcutflow=new TH1I("hcutflow","How many events still remain after cut",20,0,20);

  ShowerParticle sbjet[MAXSP],slep[MAXSP];
  Double_t HT,HT_4j,meff,meff_4j=0.,delphi,delphimin,ET_meff_4j,ET_HT_4j;
  Double_t HT_inc,MT,meff_inc,ET_HT_inc;
  Int_t njet50,nbtag50; //number of jets with pt>50

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    njet50=0;  nbtag50=0;

    progress_bar(ientry,nEvents);  //just for fun
    read_event(ientry);  //read ientry of delphes.root

    // ---------------- isolation/reconstruction/overlap removal --------------------

    //fill our histograms with reconstructed particle data
    hnelec->Fill(nelec);
    hnmuon->Fill(nmuon);
    hnlep->Fill(nelec+nmuon); 
    hnjet->Fill(njet);
    hntrack->Fill(ntrack);

    if(nelec>=1){helecpt1->Fill(selec[0].pt());}
    if(nelec>=2){helecpt2->Fill(selec[1].pt());}
    for(Int_t i=0;i<nelec;++i){
      hchargee->Fill( selec[i].charge() );
      hhademe->Fill( selec[i].hadem() );
    }

    if(nmuon>=1){hmuonpt1->Fill(smuon[0].pt());}
    if(nmuon>=2){hmuonpt2->Fill(smuon[1].pt());}
    for(Int_t i=0;i<nmuon;++i){
      hchargem->Fill( smuon[i].charge() );
    } 
   
    if(njet>=1){hjetpt1->Fill(sjet[0].pt());}
    if(njet>=2){hjetpt2->Fill(sjet[1].pt());}
    if(njet>=3){hjetpt3->Fill(sjet[2].pt());}
    if(njet>=4){hjetpt4->Fill(sjet[3].pt());}
    for(Int_t i=0;i<njet;++i){
      hmassj->Fill( sjet[i].mass() );
      hhademj->Fill( sjet[i].hadem() );
    }

    hmetpt->Fill(smissing.pt());
    hmetphi->Fill(smissing.phi());

    //baseline reconstruction
    //require jet candidates to have pt>20, |eta|<4.5
    reconstruct(sjet,njet,20.,4.5);
    //require electron candidates to have pt>20, |eta|<2.47
    reconstruct(selec,nelec,20.,2.47);
    //require muon candidates to have pt>10, |eta|<2.4
    reconstruct(smuon,nmuon,10.,2.40);

    recalculate_smissing();  //missing energy is based off reconstructed objects
    hrmetpt->Fill(smissing.pt());
    hrmetphi->Fill(smissing.phi());

    //remove jet candidates within delR=0.2 of an electron
    remove_overlap(sjet,selec,njet,nelec,0.2);
    //remove electron candidates within delR=0.4 of a jet
    remove_overlap(selec,sjet,nelec,njet,0.4);
    //remove muon candidates within delR=0.4 of a jet
    remove_overlap(smuon,sjet,nmuon,njet,0.4);

    // in this analysis an operating point is used where the btag efficiency is 0.700
    // while the light quark, c-quark and tau-lepton mistag rates are (0.001,0.200,0.077)
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.75 with light quark, c-quark and tau-lepton mistag rates of (0.017,0.250,0.125)
    nbtag=collect_bjets(sjet,sbjet,njet,20.,2.5);  //collect btagged jets in sbjet, with pt>20 and |eta|<2.5
    hnbtag->Fill(nbtag);

    //signal reconstruction
    //require jet candidates to have pt>30, |eta|<2.4
    reconstruct(sjet,njet,30.,2.4);
    //require signal muons to have pt>20, |eta|<2.40
    reconstruct(smuon,nmuon,20.,2.40);

    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon); //collect electrons and muons in slep
    hnrelep->Fill(nlep);

    // ---------------- baseline cuts ---------------------

    //require 4 jets
    if(njet<4){
      hcut->Fill(1);
      continue;
    }
    hcutflow->Fill(1);
    //requirement on leading jet pt
    if(sjet[0].pt()<90){
      hcut->Fill(2);
      continue;
    }
    hcutflow->Fill(2);
    //requirement on missing energy
    if(smissing.pt()<150){
      hcut->Fill(3);
      continue;
    }
    hcutflow->Fill(3);

    // 0 leptons?
    if(nlep==0){ 
      hcutflow->Fill(4);
      HT=getHT(sjet,njet);
      HT_4j=getHT(sjet,njet,4);
      meff=smissing.pt()+HT;
      meff_4j=smissing.pt()+HT_4j;
      ET_meff_4j=smissing.pt()/meff_4j;
      ET_HT_4j=smissing.pt()/TMath::Sqrt(HT_4j);
      
      delphimin=LONG_MAX;
      for(Int_t i=0;i<4;++i){
        delphi=sjet[i].DelPhiAbs(smissing);
        if(delphi<delphimin){
          delphimin=delphi;
        }
      } 

      hHT->Fill(HT);
      hHT_4j->Fill(HT_4j);
      hmeff->Fill(meff);
      hmeff_4j->Fill(meff_4j);
      hET_meff_4j->Fill(ET_meff_4j);
      hET_HT_4j->Fill(ET_HT_4j);
      hdelphimin->Fill(delphimin);

      if(delphimin<0.5){
        hcut->Fill(6);
        continue;
      }
      hcutflow->Fill(5);
      if(ET_meff_4j<0.2){
        hcut->Fill(7);
        continue;
      }
      hcutflow->Fill(6);
      for(Int_t i=0;i<njet;++i){
        if(sjet[i].pt()>50){
          ++njet50;
          if(sjet[i].btag()>0){
            ++nbtag50;
          }
        } 
      }
      hnjet50->Fill(njet50);
      hnbtag50->Fill(nbtag50);

      if( (smissing.pt()>200) && (nbtag>=3) && (meff_4j>1000) && (ET_HT_4j>16.) ){
        increase_SR_count(1);
      }
      if( (njet50>=4) && (nbtag50>=3) && (smissing.pt()>350) && (meff_4j>1100) ){
        increase_SR_count(2);
      }
      if( (njet50>=4) && (nbtag50>=3) && (smissing.pt()>250) && (meff_4j>1300) ){
        increase_SR_count(3);
      }  
      if( (njet>=7) && (nbtag>=3) && (smissing.pt()>200) && (meff_4j>1000) ){
        increase_SR_count(4);
      }  
      if( (njet>=7) && (nbtag>=3) && (smissing.pt()>350) && (meff_4j>1000) ){
        increase_SR_count(5);
      }  
      if( (njet>=7) && (nbtag>=3) && (smissing.pt()>250) && (meff_4j>1500) ){
        increase_SR_count(6);
      }    

    }
    else{  //or at least 1 lepton
      hcutflow->Fill(8);
      HT_inc=slep[0].pt()+getHT(sjet,njet);
      meff_inc=smissing.pt()+HT_inc;
      ET_HT_inc=smissing.pt()/TMath::Sqrt(meff_inc);
      MT=calculateMT(smissing,slep[0]);

      hHT_inc->Fill(HT_inc);
      hmeff_inc->Fill(meff_inc);
      hET_HT_inc->Fill(ET_HT_inc);
      hMT->Fill(MT);

      if( (njet>=6) && (nbtag>=3) && (smissing.pt()>175) && (MT>140) && (meff_inc>700) && (ET_HT_inc>5.) ){
        increase_SR_count(7);
      }
      if( (njet>=6) && (nbtag>=3) && (smissing.pt()>225) && (MT>140) && (meff_inc>800) && (ET_HT_inc>5.) ){
        increase_SR_count(8);
      }
      if( (njet>=6) && (nbtag>=3) && (smissing.pt()>275) && (MT>160) && (meff_inc>900) && (ET_HT_inc>5.) ){
        increase_SR_count(9);
      }  
    }
  }

  //this section outputs the efficiencies to efname.txt
  cout << "\nEfficiencies for SR_1 ... SR_" << nSR << endl;

  ofstream outFile;
  outFile.open(efname->c_str(),fstream::trunc);
  if (outFile.is_open()){
    for (int i=0;i<nSR;++i){
      outFile << efficiency[i] << ", ";
      cout << efficiency[i] << ", ";
    }
    outFile << endl;
  }
  else{
    cout << "Could not open " << efname->c_str() << " for writing.\n";
  }
  f->Write();
};

#endif
