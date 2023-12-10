#ifndef CODE_CMS2012_028_H
#define CODE_CMS2012_028_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_CMS2012_028(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnjet  = new TH1I("hnjet","Number of Jets",25,0,25);

  TH1F* hjetpt1 = new TH1F("hjetpt1","P_t of First Jet",100,0,1000);
  TH1F* hjetpt2 = new TH1F("hjetpt2","P_t of Second Jet",100,0,1000);
  TH1F* hjetpt3 = new TH1F("hjetpt3","P_t of Third Jet",100,0,1000);
  TH1F* hjetpt4 = new TH1F("hjetpt4","P_t of Fourth Jet",100,0,1000);
  TH1F* hmassj = new TH1F("hmassj","Jet's Mass",100,0,200);
  TH1F* hhademj = new TH1F("hhademj","Jet's Had/Em",100,0,20);

  TH1F* hmetpt = new TH1F("hmetpt","Missing Transverse Energy",100,0,1000);
  TH1F* hmetphi= new TH1F("hmetphi","Phi of Missing Transverse Energy",32,-TMath::Pi(),TMath::Pi());  

  TH1F* hal_t = new TH1F("hal_t","alpha_t paramater for HT>375 GeV",80,0.,4.);

  TH1I* hcut=new TH1I("hcut","Where the cuts occur",25,0,25);

  ShowerParticle sbjet[MAXSP];

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    progress_bar(ientry,nEvents);
    read_event(ientry);

    // ---------------- isolation/reconstruction/overlap removal --------------------

    //require electron candidates to have pt>10, no eta requirement
    reconstruct(selec,nelec,10.,5.0);
    //require muon candidates to have pt>10, no eta requirement
    reconstruct(smuon,nmuon,10.,5.0);
    //require photon cadidates to have pt>25, no eta requirement
    reconstruct(sphot,nphot,25.,5.0);

    // in this analysis an operating point is used where the btag efficiency is 0.6-0.7 (pt dependent)
    // while the light quark and c-quark and mistag rates are roughly 0.01
    //due to time contraints the Delphes Card that is used operates at a btag efficiency 
    // of 0.70 with light quark, c-quark mistag rates of (0.01,0.25)

    //fill our histograms with precut jet/met data
    hnjet->Fill(njet);

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

    // ---------------- cuts ---------------------

    if( (nelec+nmuon+nphot)>0 ){  //no leptons or photons
      hcut->Fill(1);
      continue;
    }

    //requirement on eta of leading jet
    if(TMath::Abs(sjet[0].eta())>2.5){
      hcut->Fill(2);
      continue;
    }
    
    Double_t htoet1,htoet2,htoet3;  //HT_slash over ET_slash
    ShowerParticle sht1,sht2,sht3;  //vector sums of all jets
    Double_t HT1=0.,HT2=0.,HT3=0.;  //scalar sum of all jet pt
    Double_t al_t;   //the alpha_t variable
    Int_t nbtag1=0,nbtag2=0,nbtag3=0;

    // here follows 58 signal regions
    //there is an interesting section regarding dedicated vetoes for inefficient
    //calorimeter cells to suppress background from energy mismeasurements
    // (last paragraph of page 7). This veto is not here replicated.

    // ----------   275 < HT < 325   ---------------------
    //require jet candidates to have pt>37 and |eta|<3.0
    reconstruct(sjet,njet,37.,3.0);
    if(njet<2){
      hcut->Fill(3);
      continue;
    }
    if( (sjet[0].Et()<73) || (sjet[1].Et()<73) ){
      hcut->Fill(4);
      continue;
    }
    for(Int_t i=0;i<njet;++i){
      HT1+=sjet[i].Et();
      sht1.Add(sjet[i]);
      if(sjet[i].btag()>0){
        ++nbtag1;
      }
    }
    if(alphat(sjet,njet)<0.55){hcut->Fill(5);} else{
      if(sht1.pt()/smissing.pt()>1.25){hcut->Fill(6);} else{
        if( (HT1<275) || (HT1>325) ){hcut->Fill(7);} else{
          if(njet<=3){
            if(nbtag1==0){increase_SR_count(1);}
            else if(nbtag1==1){increase_SR_count(9);}
            else if(nbtag1==2){increase_SR_count(17);}
          }
          else{  //njet >=4
            if(nbtag1==0){increase_SR_count(25);}
            else if(nbtag1==1){increase_SR_count(33);}
            else if(nbtag1==2){increase_SR_count(41);}
            else if(nbtag1==3){increase_SR_count(49);}    
            else if(nbtag1>=4){increase_SR_count(57);}    
          }
        }
      }
    }

    // ----------   325 < HT < 375   ---------------------
    //require jet candidates to have pt>43 and |eta|<3.0
    reconstruct(sjet,njet,43.,3.0);
    if(njet<2){
      hcut->Fill(8);
      continue;
    }
    if( (sjet[0].Et()<87) || (sjet[1].Et()<87) ){
      hcut->Fill(9);
      continue;
    }
    for(Int_t i=0;i<njet;++i){
      HT2+=sjet[i].Et();
      sht2.Add(sjet[i]);
      if(sjet[i].btag()>0){
        ++nbtag2;
      }
    }
    if(alphat(sjet,njet)<0.55){hcut->Fill(10);} else{
      if(sht2.pt()/smissing.pt()>1.25){hcut->Fill(11);} else{
        if( (HT2<325) || (HT2>375) ){hcut->Fill(12);} else{
          if(njet<=3){
            if(nbtag2==0){increase_SR_count(2);}
            else if(nbtag2==1){increase_SR_count(10);}
            else if(nbtag2==2){increase_SR_count(18);}
          }
          else{  //njet >=4
            if(nbtag2==0){increase_SR_count(26);}
            else if(nbtag2==1){increase_SR_count(34);}
            else if(nbtag2==2){increase_SR_count(42);}
            else if(nbtag2==3){increase_SR_count(50);}    
            else if(nbtag2>=4){increase_SR_count(58);}    
          }
        }
      }
    }
    // ----------   375+j*100 < HT < 475+j*100   ---------------------
    //require jet candidates to have pt>50 and |eta|<3.0
    reconstruct(sjet,njet,50.,3.0);
    if(njet<2){
      hcut->Fill(13);
      continue;
    }
    if( (sjet[0].Et()<100) || (sjet[1].Et()<100) ){
      hcut->Fill(14);
      continue;
    }
    for(Int_t i=0;i<njet;++i){
      HT3+=sjet[i].Et();
      sht3.Add(sjet[i]);
      if(sjet[i].btag()>0){
        ++nbtag3;
      }
    }
    if(sht3.pt()/smissing.pt()>1.25){
      hcut->Fill(16);
      continue;
    }
    al_t=alphat(sjet,njet);
    if(njet>=4){hal_t->Fill(al_t);}
    if(al_t<0.55){
      hcut->Fill(15);
      continue;
    }
    if((njet>=4)&&(nbtag3>=4)&&(HT3>375)){increase_SR_count(59);}    

    for(Int_t j=0;j<5;++j){
      if( (HT3<(375+j*100)) || (HT3>(475+j*100)) ){
        hcut->Fill(17+j);
        continue;
      }
      if(njet<=3){
        if(nbtag3==0){increase_SR_count(3+j);}
        else if(nbtag3==1){increase_SR_count(11+j);}
        else if(nbtag3==2){increase_SR_count(19+j);}
      }
      else{  //njet >=4
        if(nbtag3==0){increase_SR_count(27+j);}
        else if(nbtag3==1){increase_SR_count(35+j);}
        else if(nbtag3==2){increase_SR_count(43+j);}
        else if(nbtag3==3){increase_SR_count(51+j);}    
      }
    }

    // ----------   HT > 875   ---------------------
    //require jet candidates to have pt>43 and |eta|<3.0
    if( HT3<875 ){
      hcut->Fill(22);
      continue;
    }
    if(njet<=3){
      if(nbtag2==0){increase_SR_count(8);}
      else if(nbtag3==1){increase_SR_count(16);}
      else if(nbtag3==2){increase_SR_count(24);}
    }
    else{  //njet >=4
      if(nbtag3==0){increase_SR_count(32);}
      else if(nbtag3==1){increase_SR_count(40);}
      else if(nbtag3==2){increase_SR_count(48);}
      else if(nbtag3==3){increase_SR_count(56);}    
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
