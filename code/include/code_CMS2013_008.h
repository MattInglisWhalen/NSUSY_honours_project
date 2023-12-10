#ifndef CODE_CMS2013_008_H
#define CODE_CMS2013_008_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_CMS2013_008(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnelec = new TH1I("hnelec","Number of Electrons",10,0,10);
  TH1I* hnmuon = new TH1I("hnmuon","Number of Muons",10,0,10);
  TH1I* hnlep  = new TH1I("hnlep","Number of Leptons",8,0,8);
  TH1I* hnjet  = new TH1I("hnjet","Number of Jets",25,0,25);
  TH1I* hnbtag  = new TH1I("hnbtag","Number of B-tagged Jets",25,0,25);

  TH1F* helecpt1 = new TH1F("helecpt1","P_t of First Electron",100,0,1000);
  TH1F* helecpt2 = new TH1F("helecpt2","P_t of Second Electron",100,0,1000);
  TH1F* helecpt3 = new TH1F("helecpt3","P_t of Third Electron",100,0,1000);

  TH1F* hchargee = new TH1F("hchargee","Electron's Charge",20,-5,5);
  TH1F* hhademe = new TH1F("hhademe","Electron's Had/Em",100,0,0.25);

  TH1F* hmuonpt1 = new TH1F("hmuonpt1","P_t of First Muon",100,0,1000);
  TH1F* hmuonpt2 = new TH1F("hmuonpt2","P_t of Second Muon",100,0,1000);
  TH1F* hmuonpt3 = new TH1F("hmuonpt3","P_t of Third Muon",100,0,1000);
  TH1F* hchargem = new TH1F("hhchargem","Muon's Charge",20,-5,5);

  TH1F* hjetpt1 = new TH1F("hjetpt1","P_t of First Jet",100,0,1000);
  TH1F* hjetpt2 = new TH1F("hjetpt2","P_t of Second Jet",100,0,1000);
  TH1F* hjetpt3 = new TH1F("hjetpt3","P_t of Third Jet",100,0,1000);
  TH1F* hjetpt4 = new TH1F("hjetpt4","P_t of Fourth Jet",100,0,1000);
  TH1F* hmassj = new TH1F("hmassj","Jet's Mass",100,0,200);
  TH1F* hhademj = new TH1F("hhademj","Jet's Had/Em",100,0,20);

  TH1F* hmetpt = new TH1F("hmetpt","Missing Transverse Energy",100,0,1000);
  TH1F* hmetphi= new TH1F("hmetphi","Phi of Missing Transverse Energy",32,-TMath::Pi(),TMath::Pi());  

  TH1F* hHT = new TH1F("hHT","HT variable",100,0,2500);
  TH1F* hzmass=new TH1F("hzmass","Invariant mass of the reconstructed Z boson ",100,0,1000);

  TH1I* hcut=new TH1I("hcut","Where the cuts occur",20,0,20);

  Double_t HT;
  ShowerParticle slep[MAXSP],sbjet[MAXSP],stau[MAXSP];

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    progress_bar(ientry,nEvents);
    read_event(ientry);

    // ---------------- isolation/reconstruction/overlap removal --------------------

    //require jet candidates to have pt>30, |eta|<2.4
    reconstruct(sjet,njet,30.,2.4);
    //require electron candidates to have pt>10, |eta|<2.4
    reconstruct(selec,nelec,10.,2.4);
    //require muon candidates to have pt>10, |eta|<2.4
    reconstruct(smuon,nmuon,10.,2.4);

    //remove electron candidates within delR=0.1 of a muon
    remove_overlap(selec,smuon,nelec,nmuon,0.1);
    //if electron is within delR=0.4 of a jet, add it to jet
    remove_overlap(sjet,selec,njet,nelec,0.4,2);
    //if muon is within delR=0.4 of a jet, add it to jet
    remove_overlap(sjet,smuon,njet,nmuon,0.4,2);

    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon);  //put electrons and muons in slep
    // in this analysis an operating point is used where the btag efficiency is roughly 0.7
    // while the light quark and c-quark and mistag rates are less than 0.01
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.70 with light quark, c-quark mistag rates of (0.01,0.25)
    nbtag=collect_bjets(sjet,sbjet,njet);  //collect btagged jets in sbjet
    HT=getHT(sjet,njet);  //calculate HT

    //fill our histograms with precut jet/met data
    hnelec->Fill(nelec);
    hnmuon->Fill(nmuon);
    hnlep->Fill(nlep);
    hnjet->Fill(njet);
    hnbtag->Fill(nbtag);

    if(nelec>=1){helecpt1->Fill(selec[0].pt());}
    if(nelec>=2){helecpt2->Fill(selec[1].pt());}
    if(nelec>=3){helecpt3->Fill(selec[2].pt());}
    for(Int_t i=0;i<nelec;++i){
      hchargee->Fill( selec[i].charge() );
      hhademe->Fill( selec[i].hadem() );
    }

    if(nmuon>=1){hmuonpt1->Fill(smuon[0].pt());}
    if(nmuon>=2){hmuonpt2->Fill(smuon[1].pt());}
    if(nmuon>=3){hmuonpt3->Fill(smuon[2].pt());}
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

    hHT->Fill(HT);

    // ---------------- baseline cuts ---------------------

    if( nlep<3 ){  //trilepton analysis
      hcut->Fill(1);
      continue;
    }
    //require one lepton with pt>20 GeV
    if( slep[0].pt()<20 ){
      hcut->Fill(2);
      continue;
    }
    if(nbtag<1){
      hcut->Fill(3);
      continue;
    }
    if(njet<2){
      hcut->Fill(4);
      continue;
    }
    if(smissing.pt()<50){
      hcut->Fill(5);
      continue;
    }
    if(HT<60){
      hcut->Fill(6);
      continue;
    }

    ShowerParticle Zboson;
    for(Int_t i=0;i<nlep-1;++i){
      for(Int_t j=1;j<nlep;++j){
        if( (slep[i].type()==slep[j].type()) && (slep[i].sign()!=slep[j].sign()) ){
          if(Zboson.type()==UNINITIALIZED){
            Zboson=slep[i]+slep[j];
          }
          else{
            // There's another possible lepton combination to form the Z-boson 
          }
        }
      }
    }
    if(Zboson.type()!=UNINITIALIZED){
      hzmass->Fill(Zboson.mass());
      if(Zboson.mass()<12){
        hcut->Fill(7);
        continue;
      }
    }

    // ---------------- signal region definitions ---------------------

    if( (Zboson.type()!=UNINITIALIZED) && (Zboson.mass()>75) && (Zboson.mass()<105) ){  //On-Z
      if(HT<200){
        if(nbtag==1){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(1);
            }
            else if(smissing.pt()<200){
              increase_SR_count(2);
            }
            else{ // pt_miss>200
              increase_SR_count(3);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(4);
            }
            else if(smissing.pt()<200){
              increase_SR_count(5);
            }
            else{ // pt_miss>200
              increase_SR_count(6);
            }
          }
        }
        else if(nbtag==2){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(7);
            }
            else if(smissing.pt()<200){
              increase_SR_count(8);
            }
            else{ // pt_miss>200
              increase_SR_count(9);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(10);
            }
            else if(smissing.pt()<200){
              increase_SR_count(11);
            }
            else{ // pt_miss>200
              increase_SR_count(12);
            }
          }
        }
        else{ //nbtag>=3
          if(smissing.pt()<100){
            increase_SR_count(13);
          }
          else if(smissing.pt()<200){
            increase_SR_count(14);
          }
          else{ // pt_miss>200
            increase_SR_count(15);
          }
        }
      }
      else{  //>200
        if(nbtag==1){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(16);
            }
            else if(smissing.pt()<200){
              increase_SR_count(17);
            }
            else{ // pt_miss>200
              increase_SR_count(18);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(19);
            }
            else if(smissing.pt()<200){
              increase_SR_count(20);
            }
            else{ // pt_miss>200
              increase_SR_count(21);
            }
          }
        }
        else if(nbtag==2){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(22);
            }
            else if(smissing.pt()<200){
              increase_SR_count(23);
            }
            else{ // pt_miss>200
              increase_SR_count(24);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(25);
            }
            else if(smissing.pt()<200){
              increase_SR_count(26);
            }
            else{ // pt_miss>200
              increase_SR_count(27);
            }
          }
        }
        else{ //nbtag>=3
          if(smissing.pt()<100){
            increase_SR_count(28);
          }
          else if(smissing.pt()<200){
            increase_SR_count(29);
          }
          else{ // pt_miss>200
            increase_SR_count(30);
          }
        }
      } 
    }
    else{  //Off-Z
      if(HT<200){
        if(nbtag==1){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(31);
            }
            else if(smissing.pt()<200){
              increase_SR_count(32);
            }
            else{ // pt_miss>200
              increase_SR_count(33);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(34);
            }
            else if(smissing.pt()<200){
              increase_SR_count(35);
            }
            else{ // pt_miss>200
              increase_SR_count(36);
            }
          }
        }
        else if(nbtag==2){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(37);
            }
            else if(smissing.pt()<200){
              increase_SR_count(38);
            }
            else{ // pt_miss>200
              increase_SR_count(39);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(40);
            }
            else if(smissing.pt()<200){
              increase_SR_count(41);
            }
            else{ // pt_miss>200
              increase_SR_count(42);
            }
          }
        }
        else{ //nbtag>=3
          if(smissing.pt()<100){
            increase_SR_count(43);
          }
          else if(smissing.pt()<200){
            increase_SR_count(44);
          }
          else{ // pt_miss>200
            increase_SR_count(45);
          }
        }
      }
      else{  //>200
        if(nbtag==1){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(46);
            }
            else if(smissing.pt()<200){
              increase_SR_count(47);
            }
            else{ // pt_miss>200
              increase_SR_count(48);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(49);
            }
            else if(smissing.pt()<200){
              increase_SR_count(50);
            }
            else{ // pt_miss>200
              increase_SR_count(51);
            }
          }
        }
        else if(nbtag==2){
          if(njet<=3){
            if(smissing.pt()<100){
              increase_SR_count(52);
            }
            else if(smissing.pt()<200){
              increase_SR_count(53);
            }
            else{ // pt_miss>200
              increase_SR_count(54);
            }
          }
          else{ //njet>=4
            if(smissing.pt()<100){
              increase_SR_count(55);
            }
            else if(smissing.pt()<200){
              increase_SR_count(56);
            }
            else{ // pt_miss>200
              increase_SR_count(57);
            }
          }
        }
        else{ //nbtag>=3
          if(smissing.pt()<100){
            increase_SR_count(58);
          }
          else if(smissing.pt()<200){
            increase_SR_count(59);
          }
          else{ // pt_miss>200
            increase_SR_count(60);
          }
        }
      } 
    } //fi UNINITIALIZED
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
