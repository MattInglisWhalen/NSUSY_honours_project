#ifndef CODE_ATLAS2012_145_H
#define CODE_ATLAS2012_145_H

//ROOT header files
#include <TH1.h>

//C++ header files
#include <fstream>

//My header files
#include <delphesv3_class.h>

void delphesv3::code_ATLAS2012_145(){

  TFile* f = new TFile(rootfname->c_str(),"RECREATE");

  //Histos for plotting data
  TH1I* hnphot = new TH1I("hnphot","Number of Photons",25,0,25);
  TH1I* hnelec = new TH1I("hnelec","Number of Electrons",10,0,10);
  TH1I* hnmuon = new TH1I("hnmuon","Number of Muons",10,0,10);
  TH1I* hnlep  = new TH1I("hnlep","Number of Leptons",8,0,8);
  TH1I* hnjet  = new TH1I("hnjet","Number of Jets",25,0,25);
  TH1I* hnbtag = new TH1I("hnbtag","Number of total btags",15,0,15);
  TH1I* hnttag = new TH1I("hnttag","Number of total tau-tags",15,0,15);
  TH1I* hntrack  = new TH1I("hntrack","Number of Tracks",300,0,300);

  TH1F* hphotpt1 = new TH1F("hphotpt1","P_t of First Photon",100,0,500);
  TH1F* hphotpt2 = new TH1F("hphotpt2","P_t of Second Photon",100,0,500);
  TH1F* hhademp = new TH1F("hhademp","Photon's's Had/Em",100,0,0.25);
 
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

  TH1F* htrackpt = new TH1F("htrackpt","Track Pt",50,0,100);
  TH1F* htrackno = new TH1F("htrackno","Track Eta at Outer Edge",40,-5,5);
  TH1F* htrackfo = new TH1F("htrackfo","Track Phi at Outer Edge",20,-TMath::Pi(),TMath::Pi());
  TH1F* htrackx = new TH1F("htrackx","Track Vertex X Position",100,-300,300);
  TH1F* htracky = new TH1F("htracky","Track Vertex Y Position",100,-300,300);
  TH1F* htrackz = new TH1F("htrackz","Track Vertex Z Position",100,-1,1);
  TH1F* htrackxo = new TH1F("htrackxo","Track Vertex X Position at Outer Edge",100,-1500,1500);
  TH1F* htrackyo = new TH1F("htrackyo","Track Vertex Y Position at Outer Edge",100,-1500,1500);
  TH1F* htrackzo = new TH1F("htrackzo","Track Vertex Z Position at Outer Edge",100,-3500,3500);

  TH1F* hmetpt = new TH1F("hmetpt","Missing Transverse Energy",100,0,1000);
  TH1F* hmetphi= new TH1F("hmetphi","Phi of Missing Transverse Energy",32,-TMath::Pi(),TMath::Pi());  

  TH1F* hmeff_inc = new TH1F("hmeff_inc","meff_inc",150,0,3000);
  TH1F* hmeff_4j = new TH1F("hmeff_4j","meff_4j",150,0,3000);
  TH1F* hdelphimin= new TH1F("hdelphimin","Minimum Delphi Between A Leading Jet and the Missing Energy",32,-TMath::Pi(),TMath::Pi());  
  TH1F* hETmeff = new TH1F("hETmeff","ET_miss/meff",100,0,5);

  TH1I* hcut=new TH1I("hcut","Where the cuts occur",20,0,20);

  ShowerParticle slep[MAXSP],sbjet[MAXSP],stjet[MAXSP];

  for (Long64_t ientry=0; ientry<nEvents; ++ientry) { // loop over each entry in trees.root

    progress_bar(ientry,nEvents);  // just for fun
    read_event(ientry);  // reads ientry in delphes.root

    // ---------------- isolation/reconstruction/overlap removal --------------------

    //require jet candidates to have pt>20, |eta|<4.5
    reconstruct(sjet,njet,20.,4.5);
    //require electron candidates to have pt>20, |eta|<2.47
    reconstruct(selec,nelec,20.,2.47);
    //require muon candidates to have pt>10, |eta|<2.40
    reconstruct(smuon,nmuon,10.,2.40);

    //remove jets within delR<0.2 of electron
    remove_overlap(sjet,selec,njet,nelec,0.2);
    //remove electrons within delR<0.4 of jet
    remove_overlap(selec,sjet,nelec,njet,0.4);
    //remove muons within delR<0.4 of jet
    remove_overlap(smuon,sjet,nmuon,njet,0.4);

    //have to calculate the number of btagged jets
    nlep=collect_leptons(selec,smuon,slep,nelec,nmuon); //puts leptons in slep
    // in this analysis an operating point is used where the btag efficiency is 0.750
    // while the light quark, c-quark and tau-lepton mistag rates are (0.017,0.250,0.125)
    //This is the analysis that was used to set the efficiency/mistag values on the ATLAS Delphes Card
    nbtag=collect_bjets(sjet,sbjet,njet,30,2.5);  //puts btagged jets in sbjet
    nttag=collect_taus(sjet,stjet,njet);  //puts tau jets in stjet

    recalculate_smissing();  //recalculates missing energy direction after object reconstruction

    //fill our histograms with reconstructed particle data
    hnphot->Fill(nphot);
    hnelec->Fill(nelec);
    hnmuon->Fill(nmuon);
    hnlep->Fill(nlep); 
    hnjet->Fill(njet);
    hnbtag->Fill(nbtag);
    hnttag->Fill(nttag);
    hntrack->Fill(ntrack);

    if(nphot>=1){hphotpt1->Fill(sphot[0].pt());}
    if(nphot>=2){hphotpt2->Fill(sphot[1].pt());}
    for(Int_t i=0;i<nphot;++i){
      hhademp->Fill( sphot[i].hadem() );
    }

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

    for(Int_t i=0;i<ntrack;++i){
      htrackpt->Fill(strack[i].pt());
      htrackno->Fill(strack[i].etaout());
      htrackfo->Fill(strack[i].phiout());
      htrackx->Fill(strack[i].vx());
      htracky->Fill(strack[i].vy());
      htrackz->Fill(strack[i].vz());
      htrackxo->Fill(strack[i].vxo());
      htrackyo->Fill(strack[i].vyo());
      htrackzo->Fill(strack[i].vzo());
    }

    hmetpt->Fill(smissing.pt());
    hmetphi->Fill(smissing.phi());

    // ---------------- baseline cuts ---------------------

    //no electrons or muons passing reconstruction criteria
    if( nlep>0 ){
      hcut->Fill(1);
      continue;
    }
    //at least three btagged jets
    //the efficiency is 75% with "rejection factors" for (light quarks,c quarks,taus) = (58,4,8)
    if(nbtag<3){
      hcut->Fill(2);
      continue;
    }
    //at least 4 jets
    if(njet<4){   
      hcut->Fill(3);   
      continue;
    }
    //requirement on leading jet pt
    if(sjet[0].pt()<90){
      hcut->Fill(4);
      continue;
    }
    //requirement on missing energy
    if(smissing.pt()<200){
      hcut->Fill(5);
      continue;
    }

    Double_t meff_inc=smissing.pt(),meff_4j=smissing.pt(),etmeff;
    for(Int_t i=0;i<njet;++i){
      if(sjet[i].pt()>30){
        meff_inc=meff_inc+sjet[i].pt();
      }
    }
    for(Int_t i=0;i<4;++i){
      meff_4j=meff_4j+sjet[i].pt();
    }
    etmeff=smissing.pt()/meff_4j;
    hmeff_inc->Fill(meff_inc);
    hmeff_4j->Fill(meff_4j);
    hETmeff->Fill(etmeff);

    //requirement on the et/meff ratio
    if(etmeff<0.20){
      hcut->Fill(6);
      continue;
    }

    //requirement on deltes_phi between jet_i and missing energy for i=1 to 4
    Double_t delphimin=LONG_MAX,tempdelphi;
    for(Int_t i=0;i<4;++i){
      tempdelphi=smissing.DelPhiAbs(sjet[i]);
      if(tempdelphi<delphimin){
        delphimin=tempdelphi;
      }
    }
    hdelphimin->Fill(delphimin);
    if(delphimin<0.4){
      hcut->Fill(7);
      continue;
    }

    //number of jets with pt>50
    Int_t NJ=0;
    for(Int_t i=0;i<njet;++i){
      if( (sjet[i].pt()>50)&&( TMath::Abs(sjet[i].eta()) < 2.8 ) ){
        ++NJ;
      }
    }
    if(NJ<4){
      hcut->Fill(8);
      continue;
    }

    // ------------- cuts by signal region ---------------

    Int_t NJB=0;
    for(Int_t i=0;i<nbtag;++i){
      if( sbjet[i].pt()>50 ){
        ++NJB;
      }
    }

    if( (NJ>=4) && (NJB>=3) && (meff_4j>=900) ){
      increase_SR_count(1);  //increases SR count which also updates the efficiency
    }
    if( (NJ>=4) && (NJB>=3) && (meff_4j>=1100) ){
      increase_SR_count(2);
    }
    if( (NJ>=4) && (NJB>=3) && (meff_4j>=1300) ){
      increase_SR_count(3);
    }
    if( (NJ>=6) && (meff_inc>=1100) ){
      increase_SR_count(4);
    }
    if( (NJ>=6) && (meff_inc>=1300) ){
      increase_SR_count(5);
    }
    if( (NJ>=6) && (meff_inc>=1500) ){
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
