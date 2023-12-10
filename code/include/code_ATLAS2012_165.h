#ifndef CODE_ATLAS2012_165_H
#define CODE_ATLAS2012_165_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_ATLAS2012_165(){

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

  TH1F* hmeff = new TH1F("hmeff","meff",150,0,3000);
  TH1F* hETmeff = new TH1F("hETmeff","ET_miss/meff",100,0,5);
  TH1F* hMCT12 = new TH1F("hMCT12","MCT(bjet1,bjet2)",100,0,1000);
  TH1F* hHT2 = new TH1F("hHT2","HT without leading 2 jets",100,0,1000);
  TH1F* hHT3 = new TH1F("hHT3","HT without leading 3 jets",100,0,1000);

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

    //require jet candidates to have pt>20, |eta|<4.9
    reconstruct(sjet,njet,20.,4.9);
    //require electron candidates to have pt>10, |eta|<2.47
    reconstruct(selec,nelec,10.,2.47);
    //require muon candidates to have pt>10, |eta|<2.40
    reconstruct(smuon,nmuon,10.,2.40);

    //remove jets within delR<0.2 of electron
    remove_overlap(sjet,selec,njet,nelec,0.2);
    //remove electrons within delR<0.4 of jet
    remove_overlap(selec,sjet,nelec,njet,0.4);
    //remove muons within delR<0.4 of jet
    remove_overlap(smuon,sjet,nmuon,njet,0.4);

    recalculate_smissing();

    //require jet candidates to have pt>20, |eta|<2.8
    reconstruct(sjet,njet,20.,2.8);

    // in this analysis an operating point is used where the btag efficiency is 0.60
    // while the light quark, c-quark and tau-lepton mistag rates are (0.002,0.125,0.043)
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.75 with light quark, c-quark and tau-lepton mistag rates of (0.017,0.250,0.125)
    nbtag=collect_bjets(sjet,sbjet,njet,20.,2.5);  //collect btagged jets in sbjet, with pt>20 and |eta|<2.5
    hnbtag->Fill(nbtag);

    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon); //collect electrons and muons in slep
    hnrelep->Fill(nlep);

    // ---------------- baseline cuts ---------------------

    //no electrons or muons passing reconstruction criteria
    if( nlep>0 ){
      hcut->Fill(1);
      continue;
    }
    hcutflow->Fill(1);

    //require 4 jets
    if(nbtag<2){
      hcut->Fill(2);
      continue;
    }
    hcutflow->Fill(2);

    Double_t meff=smissing.pt()+getHT(sjet,njet,3);
    Double_t etmeff=smissing.pt()/meff;
    hmeff->Fill(meff);
    hETmeff->Fill(etmeff);

    if(etmeff<0.25){
      hcut->Fill(5);
      continue;
    }

    Double_t delphiJ1=sjet[0].DelPhiAbs(smissing);
    Double_t mindelphi2=LONG_MAX,tempdelphi;
    for(Int_t i=0;i<2;++i){
      tempdelphi=sjet[i].DelPhiAbs(smissing);
      if( tempdelphi < mindelphi2 ){
        mindelphi2=tempdelphi;
      }
    }
    Double_t mindelphi3=LONG_MAX;
    for(Int_t i=0;i<3;++i){
      tempdelphi=sjet[i].DelPhiAbs(smissing);
      if( tempdelphi < mindelphi3 ){
        mindelphi3=tempdelphi;
      }
    }

    //calculating Ht2 (Ht without 2 leading jets)
    Double_t Ht2=getHT(sjet,njet,njet,2);
    hHT2->Fill(Ht2);
    //calculating Ht3 (Ht without 3 leading jets)
    Double_t Ht3=getHT(sjet,njet,njet,3);
    hHT3->Fill(Ht3);

    Double_t mCT=calculateMCT(sjet[0],sjet[1]);
    hMCT12->Fill(mCT);

    if( (smissing.pt()>150.)&&(sjet[0].pt()>130.)&&(sjet[1].pt()>50.)&&(sjet[2].pt()<50.)
      &&(mindelphi2>0.4)&&(sjet[0].btag()>0.)&&(sjet[1].btag()>0.)&&(mCT>150.) ){
      increase_SR_count(1);
    }
    if( (smissing.pt()>150.)&&(sjet[0].pt()>130.)&&(sjet[1].pt()>50.)&&(sjet[2].pt()<50.)
      &&(mindelphi2>0.4)&&(sjet[0].btag()>0.)&&(sjet[1].btag()>0)&&(mCT>200.) ){
      increase_SR_count(2);
    }
    if( (smissing.pt()>150.)&&(sjet[0].pt()>130.)&&(sjet[1].pt()>50.)&&(sjet[2].pt()<50.)
      &&(mindelphi2>0.4)&&(sjet[0].btag()>0.)&&(sjet[1].btag()>0)&&(mCT>250.) ){
      increase_SR_count(3);
    }
    if( (smissing.pt()>150.)&&(sjet[0].pt()>130.)&&(sjet[1].pt()>50.)&&(sjet[2].pt()<50.)
      &&(mindelphi2>0.4)&&(sjet[0].btag()>0.)&&(sjet[1].btag()>0)&&(mCT>300.) ){
      increase_SR_count(4);
    }
    if( (smissing.pt()>200.)&&(sjet[0].pt()>60.)&&(sjet[1].pt()>60.)&&(sjet[2].pt()<50.)
      &&(mindelphi2>0.4)&&(sjet[0].btag()>0.)&&(sjet[1].btag()>0)&&(mCT>100.)&&(Ht2<50.) ){
      increase_SR_count(5);
    }
    if( (smissing.pt()>150.)&&(sjet[0].pt()>130.)&&(sjet[1].pt()>30.)&&(sjet[1].pt()<110.)&&(sjet[2].pt()>30.)
      &&(delphiJ1>2.5)&&(mindelphi3>0.4)&&(sjet[0].btag()<0.1)&&(sjet[1].btag()>0.)&&(sjet[2].btag()>0.)&&(Ht3<50.) ){
      increase_SR_count(6);
    }
    if( (smissing.pt()>250.)&&(sjet[0].pt()>150.)&&(sjet[1].pt()>30.)&&(sjet[1].pt()<110.)&&(sjet[2].pt()>30.)
      &&(delphiJ1>2.5)&&(mindelphi3>0.4)&&(sjet[0].btag()<0.1)&&(sjet[1].btag()>0.)&&(sjet[2].btag()>0.)&&(Ht3<50.) ){
      increase_SR_count(7);
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
