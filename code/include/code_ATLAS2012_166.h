#ifndef CODE_ATLAS2012_166_H
#define CODE_ATLAS2012_166_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

#ifndef WMASS
#define WMASS 80.4
#endif

void delphesv3::code_ATLAS2012_166(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnphot = new TH1I("hnphot","Number of Photons",10,0,10);
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
   
  TH1F* hHT = new TH1F("hHT","HT",150,0,3000);
  TH1F* het_sqrtHT = new TH1F("het_sqrtHT","ET_sqrt(HT)",100,0,30);
  TH1F* hMT = new TH1F("hMT","MT",100,0,1000);
  TH1F* hMJJ = new TH1F("hMJJ","Invariant mass of Wboson jj pair",100,0,1000);
  TH1F* hMJJJ = new TH1F("hMJJJ","Invariant mass of Wboson jjj pair",100,0,1000);
  TH1F* hMT2_asym = new TH1F("hMT2_asym","asymmetric MT2",100,0,1000);
  TH1F* hMT2_tau = new TH1F("hMT2_tau","MT2 assuming tau decay product",100,0,1000);
  TH1I* hcut=new TH1I("hcut","Where the cuts occur",20,0,20);
  TH1I* hcutflow=new TH1I("hcutflow","How many events still remain after cut",20,0,20);

  ShowerParticle sbjet[MAXSP],snbjet[MAXSP],slep[MAXSP];
  Int_t nnonbtag;

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    read_event(ientry);  //read ientry of delphes.root

    // ---------------- isolation/reconstruction/overlap removal --------------------

    //fill our histograms with reconstructed particle data
    hnphot->Fill(nphot);
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

    //remove jets within delR<0.2 of electron
    remove_overlap(sjet,selec,njet,nelec,0.2);

    //recalculating ETmiss
    //"Based on all electron/muon candidates, jets after overlap removal, 
    // and other energy clusters not associated to such objects"
    recalculate_smissing();

    hrmetpt->Fill(smissing.pt());
    hrmetphi->Fill(smissing.phi());

    //require reconstructed leptons to have pt>10, no requirement for eta
    reconstruct(selec,nelec,10.);
    reconstruct(smuon,nmuon,10.);  

    //leptons must be further than 0.4 from the jets
    remove_overlap(selec,sjet,nelec,njet,0.4);
    remove_overlap(smuon,sjet,nmuon,njet,0.4);

    // in this analysis an operating point is used where the btag efficiency is 0.75
    // while the light quark mistag rates are less than 0.02
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.75 with light quark, c-quark and tau-lepton mistag rates of (0.017,0.250,0.125)
    nbtag=collect_bjets(sjet,sbjet,njet,20.,2.5);  //collect btagged jets in sbjet, with pt>20 and |eta|<2.5
    nnonbtag=collect_nonbjets(sjet,snbjet,njet,20.,2.5);  //collect btagged jets in sbjet, with pt>20 and |eta|<2.5
    hnbtag->Fill(nbtag);

    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon); //collect electrons and muons in slep
    hnrelep->Fill(nlep);

    // ---------------- baseline cuts ---------------------

    //single lepton analysis
    if(nlep!=1){
      hcut->Fill(1); 
      continue;
    }
    hcutflow->Fill(1);

    //require electron to have pt>25, |eta|<2.47
    if(slep[0].type()==ELECTRON){
      if( (slep[0].pt()<25)||(TMath::Abs(slep[0].eta())>2.47) ){
        hcut->Fill(2);
        continue;
      }
    }
    else if(slep[0].type()==MUON){ //require muon candidates to have pt>25, |eta|<2.40
      if( (slep[0].pt()<25)||(TMath::Abs(slep[0].eta())>2.40) ){
        hcut->Fill(3);
        continue;
      }
    }
    hcutflow->Fill(2);

    //require 4 jets
    if(njet<4){
      hcut->Fill(4);
      continue;
    }
    hcutflow->Fill(3);

    if ( (sjet[0].pt()<80) || (sjet[1].pt()<60) || (sjet[2].pt()<40) || (sjet[3].pt()<25) ) {
      hcut->Fill(5);
      continue;
    }
    hcutflow->Fill(4);
  
    if(nbtag<1){
      hcut->Fill(6);
      continue;
    }
    hcutflow->Fill(5);

    //create JJ pair (Wboson) with MT>60 with smallest delR
    ShowerParticle Wboson,tempWboson;
    Double_t mindelphi=LONG_MAX,tempdelphi;
    Int_t ind1=-1,ind2=-1;
    //cluster together 2 jets with inv mass>60, with min delphi
    for(Int_t i=0;i<njet-1;++i){
      for(Int_t j=i+1;j<njet;++j){
        tempWboson=sjet[i]+sjet[j];
        tempdelphi=sjet[i].DelPhiAbs(sjet[j]);
        if( (tempWboson.mass()>60)&&(tempdelphi<mindelphi) ){
          ind1=i;ind2=j;
          Wboson=tempWboson;
          mindelphi=tempdelphi;
        }
      }
    }
    if( (ind1<0) || (ind2<0) ){ //could not form a Wboson
      hcut->Fill(7);
      continue;
    }
    hcutflow->Fill(6);
    Double_t mJJ=Wboson.mass();
    hMJJ->Fill(mJJ);

    ShowerParticle closest;
    mindelphi=LONG_MAX;     
    //find closest jet to wboson
    for(Int_t i=0;i<njet;++i){
      if( (i!=ind1)&&(i!=ind2) ){
        tempdelphi=sjet[i].DelPhiAbs(Wboson);
        if(tempdelphi<mindelphi){
          closest=sjet[i];
          mindelphi=tempdelphi;
        }
      }
      else{  //this jet already forms parts of the Wboson
        //do nothing
      }
    }
    
    Double_t mJJJ=(closest+Wboson).mass();
    hMJJJ->Fill(mJJJ);

    if( (mJJJ<130)||(mJJJ>205) ){
      hcut->Fill(8);
      continue;
    }
    hcutflow->Fill(7);

    // ------------- cuts by signal region ---------------

    Double_t HT=getHT(sjet,njet,4);

    Double_t et_sqrtHT=smissing.pt()/TMath::Sqrt(HT);
    hHT->Fill(HT);
    het_sqrtHT->Fill(et_sqrtHT);

    Double_t delphiJ1=sjet[0].DelPhiAbs(smissing);
    Double_t delphiJ2=sjet[1].DelPhiAbs(smissing);
    Double_t MT_lep_ptmiss=calculateMT(slep[0],smissing);
    hMT->Fill(MT_lep_ptmiss);
   
    //really not sure about this one
    //Vis1: use lep
    //Invis1: use neutrino
    //Vis2: use "tau jet" (ie highest pt non btagged jet)
    //Invis2: use neutrino
    Double_t MT2_tau=calculateMT2(slep[0],snbjet[0],smissing);
    hMT2_tau->Fill(MT2_tau);

    //from ref 78 of the ATLAS 2012_166 paper
    //arxiv: 1203.4813
    //Bai et al, "Stop the Top Background of the Stop Search", 2012

    //there are subleties here to choosing the bjets, which you need to take a careful look at
    Double_t iMT2_asym,MT2_asym=LONG_MAX;  //there are multiple possibilities, of which we choose the smallest

    if(nbtag==1){
      //Vis1: use lep+bjet
      //Invis1: use neutrino (massless)
      //Vis2: use bjet
      //Invis2: use W boson
      //must assume that one of the two highest pt non btagged jets is, in fact a bjet
      for(Int_t i=0;i<=1;++i){
      //trying nbtagged jet i
        iMT2_asym=calculateMT2(slep[0]+sbjet[0],snbjet[i],smissing,0.,WMASS);
        if(iMT2_asym<MT2_asym){MT2_asym=iMT2_asym;}
        iMT2_asym=calculateMT2(slep[0]+snbjet[i],sbjet[0],smissing,0.,WMASS);
        if(iMT2_asym<MT2_asym){MT2_asym=iMT2_asym;}
      }
    }
    if(nbtag==2){  
      //Vis1: use lep+bjet
      //Invis1: use neutrino (massless)
      //Vis2: use bjet
      //Invis2: use W boson
      //here, both bjets have been explicitly tagged
      iMT2_asym=calculateMT2(slep[0]+sbjet[0],sbjet[1],smissing,0.,WMASS);
      if(iMT2_asym<MT2_asym){MT2_asym=iMT2_asym;}
      iMT2_asym=calculateMT2(slep[0]+sbjet[1],sbjet[0],smissing,0.,WMASS);
      if(iMT2_asym<MT2_asym){MT2_asym=iMT2_asym;}
    }
    else{  //0, 3, or 4 according to the paper. I just use not 1 or 2.
      //Vis1: use lep+bjet
      //Invis1: use neutrino (massless)
      //Vis2: use bjet
      //Invis2: use W boson
      //here, we ignore btagging information, and assume that the bjets are in the 3 leading pt jets
      for(Int_t i=0;i<=1;++i){
        for(Int_t j=i+1;j<=2;++j){
          //trying jets i & j
          iMT2_asym=calculateMT2(slep[0]+sjet[i],sjet[j],smissing,0.,WMASS);
          if(iMT2_asym<MT2_asym){MT2_asym=iMT2_asym;}
          iMT2_asym=calculateMT2(slep[0]+sjet[j],sjet[i],smissing,0.,WMASS);
          if(iMT2_asym<MT2_asym){MT2_asym=iMT2_asym;}
        }
      }
    }
    hMT2_asym->Fill(MT2_asym);

    //SR 1 ( SR D )
    if( (delphiJ1>0.8) && (delphiJ2>0.8) && (smissing.pt()>275.) && (et_sqrtHT>11.) && (MT_lep_ptmiss>130.) ){
      increase_SR_count(1);
    }

    //SR2 ( SR E )
    if( (delphiJ1>0.8) && (delphiJ2>0.8) && (smissing.pt()>275.) && (et_sqrtHT>11.) && (MT_lep_ptmiss>140.) ){
      increase_SR_count(2);
    }
 
    //SR3 ( SR tN1 )
    if( (delphiJ1>0.8) && (delphiJ2>0.8) && (smissing.pt()>150.) && (et_sqrtHT>8.) && (MT_lep_ptmiss>140.) && (MT_lep_ptmiss<250.) ){
      increase_SR_count(3);
    }

    //SR4 ( SR tN2 )
    if( (delphiJ2>0.8) && (smissing.pt()>200.) && (et_sqrtHT>13.) && (MT_lep_ptmiss>140.) && (MT2_asym>170.) ){
      increase_SR_count(4);
    }

    //SR5 ( SR tN3 )
    if( (delphiJ1>0.8) && (delphiJ2>0.8) && (smissing.pt()>225.) && (et_sqrtHT>11.) && (MT_lep_ptmiss>180.) && (MT2_asym>200.) && (MT2_tau>120) ){
      increase_SR_count(5);
    }

    //SR6 ( SR bC )
    if( (delphiJ1>0.8) && (delphiJ2>0.8) && (smissing.pt()>150.) && (et_sqrtHT>7.) && (MT_lep_ptmiss>120.) ){
      increase_SR_count(6);
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
