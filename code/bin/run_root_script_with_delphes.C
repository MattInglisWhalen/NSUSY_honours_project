#include <TSystem.h>
#include <TROOT.h>

#include <iostream>
#include <string>

void run_root_script_with_delphes(char* script,char* input){

  string* run_proc=new string( ".x ");
  run_proc->append(script);
  run_proc->append(".C++(\"" );
  run_proc->append(input);
  run_proc->append( "\");" );

  std::cout<<"root [1] " << run_proc->c_str() << std::endl;

  gSystem->Load("libDelphes");
  gROOT->ProcessLine(run_proc->c_str());

  delete run_proc;

};
