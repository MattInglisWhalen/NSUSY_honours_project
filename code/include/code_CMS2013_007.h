#ifndef CODE_CMS2013_007_H
#define CODE_CMS2013_007_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>
#include <vector>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_CMS2013_007(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnelec = new TH1I("hnelec","Number of Electrons",10,0,10);
  TH1I* hnmuon = new TH1I("hnmuon","Number of Muons",10,0,10);
  TH1I* hnlep1  = new TH1I("hnlep1","Number of Leptons after pt>15 cut",8,0,8);
  TH1I* hnlep2  = new TH1I("hnlep2","Number of Leptons after pt>20 cut",8,0,8);
  TH1I* hnlep3  = new TH1I("hnlep3","Number of Leptons after Overlap Removal",8,0,8);
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

  TH1F* hHT = new TH1F("hHT","HT",50,0,2500);
  TH1F* hmetpt500= new TH1F("hmetpt500","Missing Transverse Energy for HT>500",20,0,1000);
  TH1F* hmetpt750= new TH1F("hmetpt750","Missing Transverse Energy for HT>750",20,0,1000);
  TH1F* hmetpt1000= new TH1F("hmetpt1000","Missing Transverse Energy for HT>1000",20,0,1000);
  TH1F* hdelphiWl = new TH1F("hdelphiWl","Delta Phi between the W and the lepton",12,0,TMath::Pi());
  TH1F* hdelphiWl250 = new TH1F("hdelphiWl250","Delta Phi between the W and the lepton for 250<Svar<350",12,0,TMath::Pi());
  TH1F* hdelphiWl350 = new TH1F("hdelphiWl350","Delta Phi between the W and the lepton for 350<Svar<450",12,0,TMath::Pi());
  TH1F* hdelphiWl450 = new TH1F("hdelphiWl450","Delta Phi between the W and the lepton for Svar>450",12,0,TMath::Pi());
  TH1F* hMTW = new TH1F("hMTW","Transverse mass of the lepton+smissing vector",100,0,1000);
  TH1F* hSvar = new TH1F("hSvar","The lepton scale of the event",50,0,2500);

  TH1I* hcut=new TH1I("hcut","Where the cuts occur",20,0,20);
  TH1I* hcutflow=new TH1I("hcutflow","How many events still remain after cut",30,0,30);

  ShowerParticle sWboson,sbjet[MAXSP],slep[MAXSP];
  Double_t HT,delphiWl,Svar,MT;
  Int_t nlep15;  //number of leptons with pt>15

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    progress_bar(ientry,nEvents);  //just for fun
    read_event(ientry);  //read ientry from delphes.root
    hcutflow->Fill(0);

    // ---------------- isolation/reconstruction/overlap removal --------------------

    hmetpt->Fill(smissing.pt());
    hmetphi->Fill(smissing.phi());

    hnelec->Fill(nelec);
    if(nelec>=1){helecpt1->Fill(selec[0].pt());}
    if(nelec>=2){helecpt2->Fill(selec[1].pt());}
    for(Int_t i=0;i<nelec;++i){
      hchargee->Fill( selec[i].charge() );
      hhademe->Fill( selec[i].hadem() );
    }

    hnmuon->Fill(nmuon);
    if(nmuon>=1){hmuonpt1->Fill(smuon[0].pt());}
    if(nmuon>=2){hmuonpt2->Fill(smuon[1].pt());}
    for(Int_t i=0;i<nmuon;++i){
      hchargem->Fill( smuon[i].charge() );
    } 

    //require electron candidates to have pt>15, |eta|<2.5
    reconstruct(selec,nelec,15.,2.5);
    //require muon candidates to have pt>15, |eta|<2.5
    reconstruct(smuon,nmuon,15.,2.4);
    nlep15=collect_leptons(selec,smuon,slep,nelec,nmuon);
    hnlep1->Fill(nlep15); 

    //Yes I know this is a bit weird, but signal leptons must have
    // pt>20, while requiring that no further lepton exists with pt>15,
    // which is what the above reconstruction is doing

    //require signal electrons to have pt>20, |eta|<2.5
    reconstruct(selec,nelec,20.,2.5);
    //require signal muons to have pt>20, |eta|<2.5
    reconstruct(smuon,nmuon,20.,2.4);
    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon);
    hnlep2->Fill(nlep);

/*    //remove jet candidates within delR=0.2 of an electron
    remove_overlap(sjet,selec,njet,nelec,0.2);
    //remove electron candidates within delR=0.4 of a jet
    remove_overlap(selec,sjet,nelec,njet,0.4);
    //remove muon candidates within delR=0.4 of a jet
    remove_overlap(smuon,sjet,nmuon,njet,0.4);
    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon);  */ //the overlap removal is not in the paper so I'll leave it out
    hnlep3->Fill(nlep);                                     //it also decreases the efficiency of the SRs for the signal

    if( (nlep15>=2)||(nlep!=1) ){
      hcut->Fill(2);
      continue;
    }
    hcutflow->Fill(1);

    //require jet candidates to have pt>40, |eta|<2.5
    reconstruct(sjet,njet,40.,2.4);
    // in this analysis an operating point is used where the btag efficiency is 0.7 
    // while the light quark and c-quark and mistag rates are (0.03,0.2)
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.70 with light quark, c-quark mistag rates of (0.01,0.25)
    nbtag=collect_bjets(sjet,sbjet,njet);

    recalculate_smissing();
    hrmetpt->Fill(smissing.pt());
    hrmetphi->Fill(smissing.phi());

    //fill our histograms with reconstructed particle data
    hnjet->Fill(njet);
    hntrack->Fill(ntrack);
   
    if(njet>=1){hjetpt1->Fill(sjet[0].pt());}
    if(njet>=2){hjetpt2->Fill(sjet[1].pt());}
    if(njet>=3){hjetpt3->Fill(sjet[2].pt());}
    if(njet>=4){hjetpt4->Fill(sjet[3].pt());}
    for(Int_t i=0;i<njet;++i){
      hmassj->Fill( sjet[i].mass() );
      hhademj->Fill( sjet[i].hadem() );
    }

    // ---------------- baseline cuts ---------------------

    HT=getHT(sjet,njet);
    hHT->Fill(HT);
    if(HT<500){
      hcut->Fill(3);
      continue;
    }
    hcutflow->Fill(2);

    if(njet<3){
      hcut->Fill(4);
      continue;
    }
    hcutflow->Fill(3);

    // ------------ Lepton Spectrum Method --------------

    if((njet>=6)&&(nbtag>=2)){
      hcutflow->Fill(4);
        if(HT>500){
          hcutflow->Fill(5);
          hmetpt500->Fill(smissing.pt());
          if( (smissing.pt()>=250)&&(smissing.pt()<350) ){
            increase_SR_count(1);
          }
          if( (smissing.pt()>=350)&&(smissing.pt()<450) ){
            increase_SR_count(2);
          }
          if(smissing.pt()>=450){
            increase_SR_count(3);
          }
        }
        if(HT>750){
          hcutflow->Fill(6);
          hmetpt750->Fill(smissing.pt());
          if( (smissing.pt()>=250)&&(smissing.pt()<350) ){
            increase_SR_count(4);
          }
          if( (smissing.pt()>=350)&&(smissing.pt()<450) ){
            increase_SR_count(5);
          }
          if(smissing.pt()>=450){
            increase_SR_count(6);
          }
        }
        if(HT>1000){
          hcutflow->Fill(7);
          hmetpt1000->Fill(smissing.pt());
          if( (smissing.pt()>=250)&&(smissing.pt()<350) ){
            increase_SR_count(7);
          }
          if( (smissing.pt()>=350)&&(smissing.pt()<450) ){
            increase_SR_count(8);
          }
          if(smissing.pt()>=450){
            increase_SR_count(9);
          }
        }
    }
    else{
      hcut->Fill(5);
    }

    // ------------- Delta Phi Method ------------

    sWboson=slep[0]+smissing;
    MT=calculateMT(slep[0],smissing);  //this is the good one, methinks
    hMTW->Fill(MT);
    // I dislike the fact that they don't explicitly state their definition for 
    // the transverse mass that they use to calculate their S (lepton scale) variable 
    //This paper uses the following definition for the lepton scale
    Svar=TMath::Sqrt( sWboson.pt()*sWboson.pt() + MT*MT  );
    //While their first paper witht the same analysis uses this definition
    //Svar=slep[0].pt()+smissing.pt();
    //I did some preliminary testing and both definitions give the same SR efficiencies
    //Ohhh I see, they're being sneaky, both definitions reduce to the same thing, ie they are equivalent
    hSvar->Fill(Svar);
    delphiWl=sWboson.DelPhiAbs(slep[0]);
    hdelphiWl->Fill(delphiWl);
    //sometimes these distributions have a positive slope, when they should all have negative slopes
    if( (njet>=6)&&(nbtag>=3)&&(nmuon==1)&&(Svar>=250)&&(Svar<350) ){hdelphiWl250->Fill(delphiWl);}  //for testing agreement
    else if( (njet>=6)&&(nbtag>=3)&&(nmuon==1)&&(Svar>=350)&&(Svar<450) ){hdelphiWl350->Fill(delphiWl);}  //for testing agreement
    else if( (njet>=6)&&(nbtag>=3)&&(nmuon==1)&&(Svar>=450) ){hdelphiWl450->Fill(delphiWl);}  //for testing agreement
    if(delphiWl<1){
      hcut->Fill(7);
      continue;
    }
    hcutflow->Fill(8);

    if(njet<=5){
      hcutflow->Fill(9);
      if(nbtag==2){
        hcutflow->Fill(10);
        if(nmuon==1){
          hcutflow->Fill(11);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(10);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(11);
          }
          else if(Svar>=450){
            increase_SR_count(12);
          }
          else{
            hcut->Fill(8);
          }
        }
        else{  //nelec==1
          hcutflow->Fill(12);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(13);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(14);
          }
          else if(Svar>=450){
            increase_SR_count(15);
          }
          else{
            hcut->Fill(9);
          }
        }
      }
      else if(nbtag>=3){
        hcutflow->Fill(13);
        if(nmuon==1){
          hcutflow->Fill(14);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(16);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(17);
          }
          else if(Svar>=450){
            increase_SR_count(18);
          }
          else{
            hcut->Fill(10);
          }
        }
        else{  //nelec==1
          hcutflow->Fill(15);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(19);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(20);
          }
          else if(Svar>=450){
            increase_SR_count(21);
          }
          else{
            hcut->Fill(11);
          }
        }
      }
      else{
        //do nothing
      }
    }
    else{ // njet>=6 
      hcutflow->Fill(16);
      if(nbtag==2){
        hcutflow->Fill(17);
        if(nmuon==1){
          hcutflow->Fill(18);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(22);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(23);
          }
          else if(Svar>=450){
            increase_SR_count(24);
          }
          else{
            hcut->Fill(12);
          }
        }
        else{  //nelec==1
          hcutflow->Fill(19);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(25);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(26);
          }
          else if(Svar>=450){
            increase_SR_count(27);
          }
          else{
            hcut->Fill(13);
          }
        }
      }
      else if(nbtag>=3){
        hcutflow->Fill(20);
        if(nmuon==1){
          hcutflow->Fill(21);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(28);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(29);
          }
          else if(Svar>=450){
            increase_SR_count(30);
          }
          else{
            hcut->Fill(14);
          }
        }
        else{  //nelec==1
          hcutflow->Fill(22);
          if( (Svar>=250)&&(Svar<350) ){
            increase_SR_count(31);
          }
          else if( (Svar>=350)&&(Svar<450) ){
            increase_SR_count(32);
          }
          else if(Svar>=450){
            increase_SR_count(33);
          }
          else{
            hcut->Fill(15);
          }
        }
      }
      else{
        //do nothing
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
