#ifndef CODE_ATLAS2013_007_H
#define CODE_ATLAS2013_007_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_ATLAS2013_007(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnelec = new TH1I("hnelec","Number of Electrons",10,0,10);
  TH1I* hnmuon = new TH1I("hnmuon","Number of Muons",10,0,10);
  TH1I* hnlep  = new TH1I("hnlep","Number of Leptons",8,0,8);
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

  TH1F* hsign= new TH1F("hsign","Lepton sign",10,-2,2);

  TH1F* hmeff = new TH1F("hmeff","Effective Mass",100,0,1000);
  TH1F* hMT = new TH1F("hMT","Transverse Mass",100,0,1000);

  TH1I* hcut=new TH1I("hcut","Where the cuts occur",20,0,20);
  TH1I* hcutflow=new TH1I("hcutflow","How many events still remain after cut",20,0,20);

  ShowerParticle sbjet[MAXSP],slep[MAXSP];
  Double_t MT,meff;

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    progress_bar(ientry,nEvents);  //just for fun
    read_event(ientry);  //reads ientry of delphes.root

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

    //require jet candidates to have pt>20, |eta|<2.8
    reconstruct(sjet,njet,20.,2.8);
    //the interpretation is a bit weird, but I think it says that
    // b-tagged jets only need pt>20 while non b-tagged jets
    // are required to have pt>40
    // in this analysis an operating point is used where the btag efficiency is 0.7
    // while the light quark, c-quark and tau-lepton mistag rates are undefined
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.75 with light quark, c-quark and tau-lepton mistag rates of (0.017,0.250,0.125)
    nbtag=collect_bjets(sjet,sbjet,njet);  //collect btagged jets in sbtag
    reconstruct(sjet,njet,40.,2.8);

    hnbtag->Fill(nbtag);

    //require electron candidates to have pt>20, |eta|<2.47
    reconstruct(selec,nelec,20.,2.47);
    //require muon candidates to have pt>20, |eta|<2.40
    reconstruct(smuon,nmuon,20.,2.40);

    //remove jet candidates within delR=0.2 of an electron
    remove_overlap(sjet,selec,njet,nelec,0.2,1);
    //remove electron candidates within delR=0.4 of a jet
    remove_overlap(selec,sjet,nelec,njet,0.4,1);
    //remove muon candidates within delR=0.4 of a jet
    remove_overlap(smuon,sjet,nmuon,njet,0.4,1);

    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon);  //collect leptons (elec and mu) in slep

    //recalculate missing transverse momentum
    //because this analysis only uses reconstructed objects
    recalculate_smissing();
    hrmetpt->Fill(smissing.pt());
    hrmetphi->Fill(smissing.phi());

    // ---------------- baseline cuts ---------------------

    //require primary vertex with at least 5 tracks
    //can be done rudimentarily by just requiring >=5 tracks
    if (ntrack<5){
      hcut->Fill(1);
      continue;
    }

    //dilepton analysis
    if(nlep<2){
      hcut->Fill(2);
      continue;
    }
    hcutflow->Fill(2);
    hsign->Fill(slep[0].sign());
    //hsign->Fill(slep[1].sign());
    //same sign lepton analysis
    if( slep[0].sign() != slep[1].sign() ){
      hcut->Fill(3);
      continue;
    }
    hcutflow->Fill(3);

    //fairly rudimentary analysis variables
    MT=calculateMT(smissing,slep[0]);
    meff=smissing.pt()+getHT(sjet,njet)+slep[0].pt()+slep[1].pt();

    hMT->Fill(MT);
    hmeff->Fill(meff);

    // ------------- cuts by signal region ---------------

    if (njet<3){
      hcut->Fill(4);
      continue;
    }
    hcutflow->Fill(4);

    //two signal regions which are outlined in the paper are not used here
    // because they incorporate a "binned shape fit" which is not described in the paper

    if( nbtag==0 ){
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      hcutflow->Fill(5);
      if (smissing.pt()>150) {
        hcutflow->Fill(6);
        if (MT>100) {
          hcutflow->Fill(7);
          if (meff>400) {
            hcutflow->Fill(8);
            increase_SR_count(1);
          }
        }
      }
    }
    if(nbtag>=1){
      hcutflow->Fill(9);
      if (smissing.pt()>150){
        hcutflow->Fill(10);
        if(MT>100){
          hcutflow->Fill(11);
          if(meff>700){
            hcutflow->Fill(12);
            increase_SR_count(2);
          }
        }
      }
    }
    if(nbtag>=3){
      hcutflow->Fill(13);
      if(njet>=4){
        hcutflow->Fill(14);
        increase_SR_count(3);
      }
    }
    if( (nbtag>=3) && (njet>=5) && (smissing.pt()<150) && (MT<100) ){
      increase_SR_count(4);
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
