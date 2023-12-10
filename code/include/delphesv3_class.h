#ifndef DELPHESV3_CLASS_H
#define DELPHESV3_CLASS_H

#include <string>

#include "MIW_ShowerParticle.h"

#define NUM_CODES 8
#define MAXSP 300
#define MAXTRACK 1000
#define MAXPART 2000

class delphesv3{

private:

  //analysis codes via function pointers
  typedef void (delphesv3::*code_ptr)();
  code_ptr codes[NUM_CODES];

  TChain* fChain;
  ExRootTreeReader* treeReader;
  Long64_t nEvents;

  //For this particular analysis
  Int_t nSR;
  Int_t* SR_count;
  Double_t* efficiency;
  string *ansysname,*rootfname,*efname;
  string* delphes_codelist[NUM_CODES];

  TClonesArray* b_track;
  TClonesArray* b_phot;
  TClonesArray* b_elec;
  TClonesArray* b_muon;
  TClonesArray* b_tau;
  TClonesArray* b_neut;
  TClonesArray* b_jet;
  TClonesArray* b_met;
  TClonesArray* b_part;

  //The event's shower particles
  ShowerParticle   spart[MAXPART];
  ShowerParticle     sphot[MAXSP];
  ShowerParticle     selec[MAXSP];
  ShowerParticle     smuon[MAXSP];
  ShowerParticle      stau[MAXSP];
  ShowerParticle     sneut[MAXSP];
  ShowerParticle      sjet[MAXSP];
  ShowerParticle strack[MAXTRACK];
  ShowerParticle         smissing;

  //Number of particles in each event
  Int_t  npart; //number of particles
  Int_t  nphot; //number of photons
  Int_t  nelec; //number of electrons
  Int_t  nmuon; //number of muons
  Int_t   ntau; //number of tauons
  Int_t  nneut; //number of neutrinos
  Int_t   nlep; //number of leptons (electrons+muons)
  Int_t   njet; //number of jets
  Int_t ntrack; //number of tracks
  Int_t  nbtag; //number of jets that are bottom-tagged
  Int_t  nttag; //number of jets that are tau-tagged
  Int_t   nmet; //number of "missing energy" particles (should only be 1)

  //member functions
  void read_event(Long64_t);
  Int_t get_code_index();
  void recalculate_smissing();
  void increase_SR_count(Int_t);
  void trim_delphes();
  void look_at_HEP_events();

  //analysis codes
  void code_ATLAS2012_145();
  void code_ATLAS2012_165();
  void code_ATLAS2012_166();
  void code_ATLAS2013_007();
  void code_ATLAS2013_061();
  void code_CMS2012_028();
  void code_CMS2013_007();
  void code_CMS2013_008();

 public:

  //constructor/destructors
  delphesv3(char*,TTree* = NULL);
  ~delphesv3();

  //Running analyses
  void SetAnalysis(char*);
  void RunAnalysis();
  void get_SR_eff(Int_t);

};

#endif
