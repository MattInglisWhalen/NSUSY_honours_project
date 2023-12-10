#ifndef DELPHESV3_FUNCS_H
#define DELPHESV3_FUNCS_H

#include <TChain.h>
#include <TClonesArray.h>

#include <string>
#include <iostream>

#include <classes/DelphesClasses.h>
#include <ExRootAnalysis/ExRootTreeReader.h>

#include <delphesv3_class.h>
#include <MIW_analysis_lib.h>
#include <code_ATLAS2012_145.h>
#include <code_ATLAS2012_165.h>
#include <code_ATLAS2012_166.h>
#include <code_ATLAS2013_007.h>
#include <code_ATLAS2013_061.h>
#include <code_CMS2012_028.h>
#include <code_CMS2013_007.h>
#include <code_CMS2013_008.h>
#include <look_at_HEP_events.h>

#define DEBUG 0

using namespace std;

delphesv3::delphesv3(char* ansystag,TTree* tree) : fChain(0) {

  ansysname=new string();
  rootfname=new string();
  efname=new string();
  SetAnalysis(ansystag);

  fChain=new TChain("Delphes");
  if( ansysname->find("ATLAS") != ansysname->npos ){
    fChain->Add("tag_1_delphes_events_ATLAS.root");
  }
  else if( ansysname->find("CMS") != ansysname->npos ){
    fChain->Add("tag_1_delphes_events_CMS.root");
  }
  else{
    fChain->Add("tag_1_delphes_events.root");
  }

  treeReader=new ExRootTreeReader(fChain);
  nEvents= treeReader->GetEntries();

  b_track = treeReader->UseBranch("Track");
  b_phot  = treeReader->UseBranch("Photon");
  b_elec  = treeReader->UseBranch("Electron");
  b_muon  = treeReader->UseBranch("Muon");
  b_jet   = treeReader->UseBranch("Jet");
  b_met   = treeReader->UseBranch("MissingET");
  b_part  = treeReader->UseBranch("Particle");
 
  //list of analysis code names;
  delphes_codelist[0]=new string("ATLAS2012_145");
  delphes_codelist[1]=new string("ATLAS2012_165");
  delphes_codelist[2]=new string("ATLAS2012_166");
  delphes_codelist[3]=new string("ATLAS2013_007");
  delphes_codelist[4]=new string("ATLAS2013_061");
  delphes_codelist[5]=new string("CMS2012_028");
  delphes_codelist[6]=new string("CMS2013_007");
  delphes_codelist[7]=new string("CMS2013_008");

  //initialization of function pointer list
  codes[0]=&delphesv3::code_ATLAS2012_145;
  codes[1]=&delphesv3::code_ATLAS2012_165;
  codes[2]=&delphesv3::code_ATLAS2012_166;
  codes[3]=&delphesv3::code_ATLAS2013_007;
  codes[4]=&delphesv3::code_ATLAS2013_061;
  codes[5]=&delphesv3::code_CMS2012_028;
  codes[6]=&delphesv3::code_CMS2013_007;
  codes[7]=&delphesv3::code_CMS2013_008;

};

delphesv3::~delphesv3(){
  //destructor
};

void delphesv3::read_event(Long64_t entry){

  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);

  if (DEBUG) { if (entry>0) exit(1); }

  if(b_phot&&b_elec&&b_muon&&b_jet&&b_track&&b_met){      // the ROOT file was made by Delphes
    nphot=b_phot->GetEntriesFast();
    nelec=b_elec->GetEntriesFast();
    nmuon=b_muon->GetEntriesFast();
    nmuon=b_tau->GetEntriesFast();
    nneut=b_neut->GetEntriesFast();
    njet=b_jet->GetEntriesFast();
    ntrack=b_track->GetEntriesFast(); 
    nmet=b_track->GetEntriesFast();

    for(Int_t i=0;i<nphot;++i){
      sphot[i].SetType(PHOTON);
      sphot[i].SetPtEtaPhiMass( ((Photon*)b_phot->At(i))->PT,
                                ((Photon*)b_phot->At(i))->Eta,
                                ((Photon*)b_phot->At(i))->Phi, 0. );
      sphot[i].SetChargeHademBtagTtag(0.,((Photon*)b_phot->At(i))->EhadOverEem,0.,0. );
    }
  
    for(Int_t i=0;i<nelec;++i){
      selec[i].SetType(ELECTRON);
      selec[i].SetPtEtaPhiMass( ((Electron*)b_elec->At(i))->PT,
                                ((Electron*)b_elec->At(i))->Eta,
                                ((Electron*)b_elec->At(i))->Phi,MASS_ELECTRON );
      selec[i].SetChargeHademBtagTtag( ((Electron*)b_elec->At(i))->Charge,
                                       ((Electron*)b_elec->At(i))->EhadOverEem,0.,0.);
    }

    for(Int_t i=0;i<nmuon;++i){
      smuon[i].SetType(MUON);
      smuon[i].SetPtEtaPhiMass( ((Muon*)b_muon->At(i))->PT,
                                ((Muon*)b_muon->At(i))->Eta,
                                ((Muon*)b_muon->At(i))->Phi,MASS_MUON);
      smuon[i].SetChargeHademBtagTtag(((Muon*)b_muon->At(i))->Charge,0.,0.,0.);
    }
    nlep=0;

    for(Int_t i=0;i<njet;++i){
      sjet[i].SetType(JET);
      sjet[i].SetPtEtaPhiMass(((Jet*)b_jet->At(i))->PT,
                              ((Jet*)b_jet->At(i))->Eta,
                              ((Jet*)b_jet->At(i))->Phi,
                              ((Jet*)b_jet->At(i))->Mass ) ;
      sjet[i].SetChargeHademBtagTtag(((Jet*)b_jet->At(i))->Charge,
                                     ((Jet*)b_jet->At(i))->EhadOverEem,
                                     ((Jet*)b_jet->At(i))->BTag,
                                     ((Jet*)b_jet->At(i))->TauTag );
    }
    nbtag=0;
    nttag=0;

    for(Int_t i=0;i<ntrack;++i){
      strack[i].SetType(TRACK);
      strack[i].SetPtEtaPhiMass( ((Track*)b_track->At(i))->PT,
                                                    ((Track*)b_track->At(i))->Eta,
                                                    ((Track*)b_track->At(i))->Phi,MASS_ELECTRON);
      strack[i].SetChargeHademBtagTtag(((Track*)b_track->At(i))->Charge,0.,0.,0.);
      strack[i].SetEtaOPhiOVxVyVzVxoVyoVzo(((Track*)b_track->At(i))->EtaOuter,
                                     ((Track*)b_track->At(i))->PhiOuter,
                                     ((Track*)b_track->At(i))->X,
                                     ((Track*)b_track->At(i))->Y,
                                     ((Track*)b_track->At(i))->Z,
                                     ((Track*)b_track->At(i))->XOuter,
                                     ((Track*)b_track->At(i))->YOuter,
                                     ((Track*)b_track->At(i))->ZOuter );
    }
 
    if( nmet>0 ){  //there's always only one "Missing Energy" particle
      smissing.SetPtEtaPhiMass( ((MissingET*)b_met->At(0))->MET,0.,
                               ((MissingET*)b_met->At(0))->Phi,0.  );
    }

  } else {                              // HEP file type converted to ROOT tree 
					// using ./DelphesSTDHEP examples/converter_card.tcl output.root input.hep
    nphot=nelec=nmuon=ntau=nneut=njet=ntrack=0;
    smissing.clear();
    npart=b_part->GetEntriesFast();
    if(npart>MAXPART){
      cout << "Maximum number of particles in a single event exceeded, exiting...\n";
      exit(-1);
    }

    for(Int_t i=0;i<npart;++i){   
      if( (nphot>MAXSP) || (nelec>MAXSP) || (nmuon>MAXSP) || (ntau>MAXSP) || (nneut>MAXSP) ){
        cout << "Maximum number of shower particles of a single type in a single event exceeded, exiting...\n";
        exit(-1);
      }  
      if(ntrack>MAXTRACK){
        cout << "Maximum number of tracks in a single event exceeded, exiting...\n";
        exit(-1);
      }
      if(  ((GenParticle*)b_part->At(i))->Status != 1 ){     // initial or intermediate state particle
        continue;
      }
      // else { //Final State Particle (or, this should be the way)
      spart[i].SetType( ((GenParticle*)b_part->At(i))->PID );
      spart[i].SetPtEtaPhiMass( ((GenParticle*)b_part->At(i))->PT,
                                ((GenParticle*)b_part->At(i))->Eta,
                                ((GenParticle*)b_part->At(i))->Phi, 
                                ((GenParticle*)b_part->At(i))->Mass  );
      spart[i].SetChargeHademBtagTtag(((GenParticle*)b_part->At(i))->Charge,0.,0.,0. );

      if( TMath::Abs( spart[i].type() ) == 1 ){ // down quark
        cout << "i=" << i << ": down quark in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 2 ){ // up quark
        cout << "i=" << i << ": up quark in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 3 ){ // strange quark
        cout << "i=" << i << ": strange quark in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 4 ){ // charm quark
        cout << "i=" << i << ": charm quark in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 5 ){ // bottom quark
        cout << "i=" << i << ": bottom quark in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 6 ){ // top quark
        cout << "i=" << i << ": top quark in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 11 ){ // electron
        selec[nelec++]=spart[i];
        selec[nelec-1].SetType(ELECTRON);
        if(DEBUG){ cout << i << " e " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 12 ){ // electron neutrino
        sneut[nneut++]=spart[i];
        sneut[nneut-1].SetType(INVISIBLE);
        if(DEBUG){ cout << i << " nu e " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 13 ){ // muon
        smuon[nmuon++]=spart[i];
        smuon[nmuon-1].SetType(MUON);
        if(DEBUG){ cout << i << " mu " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 14 ){ // muon neutrino
        sneut[nneut++]=spart[i];
        sneut[nneut-1].SetType(INVISIBLE);
        if(DEBUG){ cout << i << " nu mu " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 15 ){ // tauon
        stau[ntau++]=spart[i];
        stau[ntau-1].SetType(TAUON);
        if(DEBUG){ cout << i << " tau " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 16 ){ // tau neutrino
        sneut[nneut++]=spart[i];
        sneut[nneut-1].SetType(INVISIBLE);
        if(DEBUG){ cout << i << " nu tau " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 21 ){ // gluon
        cout << "i=" << i << ": gluon in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 22 ){ // photon
        sphot[nphot++]=spart[i];
        sphot[nphot-1].SetType(PHOTON);
        if(DEBUG){ cout << i << " gamma " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 23 ){ // Z boson
        cout << "i=" << i << ": Z boson in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 24 ){ // W boson
        cout << "i=" << i << ": W boson in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 25 ){ // Higgs boson
        cout << "i=" << i << ": Higgs boson in the final state.\n" << endl;
        exit(-1);
      }
      else if( TMath::Abs( spart[i].type() ) == 111 ){ // neutral pion
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " pi^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 113 ){ // neutral rho meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " rho^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 130 ){ // K long meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " K_L^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 211 ){ // charged pion
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " pi^+- " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 213 ){ // charged rho meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " rho^+- " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 221 ){ // eta meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " eta " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 223 ){ // omega meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " omega " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 310 ){ // K short meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " K_S^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 311 ){ // neutral kaon
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " K^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 313 ){ // neutral kaon (excited)
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " K^*0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 321 ){ // charged kaon
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " K^+- " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 321 ){ // charged kaon (excited)
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " K^*+- " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 331 ){ // eta prime meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " eta' " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 333 ){ // phi meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " phi " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 411 ){ // charged D meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " D_+- " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 421 ){ // neutral D meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " D_0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 431 ){ // charged strange D meson
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " D_0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 2112 ){ // neutron
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " n " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 2212 ){ // proton
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " p " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 2212 ){ // neutral Lambda baryon
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " Lambda^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 2212 ){ // negative Sigma baryon
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " Sigma^- " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 2212 ){ // neutral Sigma baryon
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " Sigma^0 " << endl;}
      }
      else if( TMath::Abs( spart[i].type() ) == 2212 ){ // positive Sigma baryon
        strack[ntrack++]=spart[i];
        if(DEBUG){ cout << i << " Sigma^+ " << endl;}
      }
      else{
        cout << "\nUnclassified particle " << spart[i].type() << "  -- please extend particle dictionary.\n";
        // Can extend dictionary via http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
        exit(-1);
      }
    }  // end of loop over particles in event

    // cluster the hadronic tracks into jets
    cluster_antikt(sjet,njet,strack,ntrack);

    // sort shower particles by pt, highest first in list
    sort_by_pt(sphot,nphot);
    sort_by_pt(selec,nelec);
    sort_by_pt(smuon,nmuon);
    sort_by_pt(stau,ntau);
    sort_by_pt(sneut,nneut);
    sort_by_pt(sjet,njet);
  
    // recalculate the missing energy vector
    recalculate_smissing();

  }  // end of if statement determining proper way to read the tree
};

Int_t delphesv3::get_code_index(){
  for(Int_t i=0;i<NUM_CODES;++i){
    if( strcmp(ansysname->c_str(),delphes_codelist[i]->c_str()) == 0 ){
      return i;
    }
  }
  return -1;
};


//public functions
void delphesv3::SetAnalysis(char* ansystag){
  ansysname->clear(); rootfname->clear(); efname->clear();
  ansysname->append(ansystag); 
  if ( ansysname->find(".root") != ansysname->npos ){
    rootfname->append("trimmed_delphes.root");
    return;
  }
  else{
    rootfname->append(ansystag);
    rootfname->append("_analysis.root");
    efname->append(ansystag);
    efname->append("_eff.txt");
  }
  nSR=get_nSR(ansystag);
  if(nSR<0){
    cout << *ansysname << " does not exist in MIW_search_data.h, please add it before continuing.\n";
    exit(-1);
  }

  delete[] SR_count;
  SR_count=new Int_t[nSR];
  memset(SR_count,0,sizeof(Int_t)*nSR);

  delete[] efficiency;
  efficiency=new Double_t[nSR];
  memset(efficiency,0,sizeof(Double_t)*nSR);
}

void delphesv3::RunAnalysis(){
  Int_t code_index=get_code_index();
  if(code_index>=0){
    ( this->*codes[code_index])(); //run the analysis code
  }
  else{
    cout << "Analysis code " << (*ansysname) << " does not exist.\n";
    exit(-1);
  }
};

void delphesv3::recalculate_smissing(){
  smissing.clear();
  for(Int_t i=0;i<nphot;++i){
    smissing.NegVecSum(sphot[i]);
  }
  for(Int_t i=0;i<nelec;++i){
    smissing.NegVecSum(selec[i]);
  }
  for(Int_t i=0;i<nmuon;++i){
    smissing.NegVecSum(smuon[i]);
  }
  for(Int_t i=0;i<ntau;++i){
    smissing.NegVecSum(stau[i]);
  }
  for(Int_t i=0;i<nneut;++i){        //This should be taken out for detector sims
    smissing.NegVecSum(sneut[i]);
  }
  for(Int_t i=0;i<njet;++i){
    smissing.NegVecSum(sjet[i]);
  }
  for(Int_t i=0;i<ntrack;++i){       // For testing purposes only
    //smissing.NegVecSum(strack[i]);
  }
  smissing.SetType(MISSING);
};

void delphesv3::increase_SR_count(Int_t SR){
  SR_count[SR-1]=SR_count[SR-1]+1;
  efficiency[SR-1]=((Double_t)SR_count[SR-1])/((Double_t)nEvents);
};

void delphesv3::trim_delphes(){

  TFile* trimmedFile=new TFile(rootfname->c_str(),"RECREATE");
  TTree* trimmedTree=new TTree("Delphes","Trimmed Analysis");

  TBranch* flag=NULL;

  flag=trimmedTree->Branch("Track",b_track);
  if(!flag){std::cout << "Track - Branch failure.\n";}

  flag=trimmedTree->Branch("Jet",b_jet);
  if(!flag){std::cout << "Jet - Branch failure.\n";}

  flag=trimmedTree->Branch("Photon",b_phot);
  if(!flag){std::cout << "Photon - Branch failure.\n";}

  flag=trimmedTree->Branch("Electron",b_elec);
  if(!flag){std::cout << "Electron - Branch failure.\n";}

  flag=trimmedTree->Branch("Muon",b_muon);
  if(!flag){std::cout << "Muon - Branch failure.\n";}

  flag=trimmedTree->Branch("MissingET",b_met);
  if(!flag){std::cout << "MissingET - Branch failure.\n";}

  flag=trimmedTree->Branch("Particle",b_part);
  if(!flag){std::cout << "Particle - Branch failure.\n";}

  for(Long64_t iEntry=0;iEntry<nEvents;++iEntry){
    read_event(iEntry);
    trimmedTree->Fill();
  }
  trimmedFile->Write();
  trimmedFile->Close();

};

#endif
