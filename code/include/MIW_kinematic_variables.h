#ifndef MIW_KINEMATIC_VARIABLES
#define MIW_KINEMATIC_VARIABLES

//root files necessary for function minimization, for MT2 in particular
#include <TFitter.h>
#include <TVirtualFitter.h>
#include <TMinuit.h>

//C++ header files
#include <vector>

//MIW header files
#include "MIW_useful_functions.h"
#include "MIW_ShowerParticle.h"

Double_t calculateMT(ShowerParticle s1, ShowerParticle s2){
  Double_t MT_squared=s1.mass()*s1.mass()+s2.mass()*s2.mass()+2*(s1.Et()*s2.Et()-s1.TransverseInnerProduct(s2));
  return TMath::Sqrt(MT_squared);
};

Double_t calculateMCT(ShowerParticle s1, ShowerParticle s2){
  Double_t MCT_squared=s1.mass()*s1.mass()+s2.mass()*s2.mass()+2*(s1.Et()*s2.Et()+s1.TransverseInnerProduct(s2));
  return TMath::Sqrt(MCT_squared);
};

namespace MIW_MT2{
  ShowerParticle vis_s1,vis_s2,missing_s;
  Double_t invis_m1,invis_m2;
}

//See H. Cheng and Z. Han. "Minimal Kinematic Constraints and $m_{T_{2}}$"
//here, we assume that the invisible momenta have zero mass
void minuit_MT2(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/  ){
  ShowerParticle q1,q2; //the invisible momenta
  Double_t MTi,MTii,px1,px2,py1,py2;
  px1=p[0]; 
  py1=p[1];
  px2=MIW_MT2::missing_s.px()-px1;
  py2=MIW_MT2::missing_s.py()-py1;
  q1.SetPtEtaPhiMass(TMath::Sqrt(px1*px1+py1*py1),0,TMath::ATan2(py1,px1),MIW_MT2::invis_m1);
  q2.SetPtEtaPhiMass(TMath::Sqrt(px2*px2+py2*py2),0,TMath::ATan2(py2,px2),MIW_MT2::invis_m2);
  MTi=calculateMT(q1,MIW_MT2::vis_s1);
  MTii=calculateMT(q2,MIW_MT2::vis_s2);
  if ( MTi>=MTii ){
    fval=MTi;
  }
  else if( MTi<MTii ){
    fval=MTii;
  }
  else{  //something very wrong happened, make sure minuit doesn't keep it
    fval=LONG_MAX;  //usually this is because either MTi or MTii is nan, don't ask me why
  }
};

Double_t calculateMT2(ShowerParticle s1, ShowerParticle s2, ShowerParticle missing_p, Double_t invis_mass1=0, Double_t invis_mass2=0){

  //the namespace is the only way to pass variables to
  //the minuit function (I don't like global variables,
  //and classes seem like too much structure)
  MIW_MT2::vis_s1=s1;
  MIW_MT2::vis_s2=s2;
  MIW_MT2::missing_s=missing_p;
  MIW_MT2::invis_m1=invis_mass1;
  MIW_MT2::invis_m2=invis_mass2;

  //Minimize the minuit_Mt2 function with Minuit and thus find the best estimator for invis_px and invis_py
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuit = TVirtualFitter::Fitter(0,2);
  minuit->SetParameter(0,"px of invisible particle 1",0.0, 50.0, 0.0, 0.0);
  minuit->SetParameter(1,"py of invisible particle 2",0.0, 50.0, 0.0, 0.0);
  minuit->SetFCN(minuit_MT2);

  // minimize myFcn
  Double_t arglist[100];
  arglist[0]=-1; //no printing!!!
  minuit->ExecuteCommand("SET PRINT",arglist,1);

  arglist[0]=5000;  //max number of function calls
  arglist[1]=0.001; //tolerance
  minuit->ExecuteCommand("SET NOWARNING",arglist,0);
  minuit->ExecuteCommand("SIMPLEX",arglist,0);
  minuit->ExecuteCommand("MIGRAD",arglist,0);
  minuit->ExecuteCommand("MINOS",arglist,0);

  Double_t optimal_invis_x1, optimal_invis_y1;
  Double_t optimal_invis_x2, optimal_invis_y2;

  optimal_invis_x1=minuit->GetParameter(0);
  optimal_invis_y1=minuit->GetParameter(1);
  optimal_invis_x2=missing_p.px()-optimal_invis_x1;
  optimal_invis_y2=missing_p.py()-optimal_invis_y1;

  ShowerParticle optimal_invis_part1,optimal_invis_part2;

  optimal_invis_part1.SetType(INVISIBLE);
  optimal_invis_part2.SetType(INVISIBLE);
  optimal_invis_part1.SetPtEtaPhiMass(TMath::Sqrt(optimal_invis_x1*optimal_invis_x1+optimal_invis_y1*optimal_invis_y1),0,TMath::ATan2(optimal_invis_y1,optimal_invis_x1),invis_mass1);
  optimal_invis_part2.SetPtEtaPhiMass(TMath::Sqrt(optimal_invis_x2*optimal_invis_x2+optimal_invis_y2*optimal_invis_y2),0,TMath::ATan2(optimal_invis_y2,optimal_invis_x2),invis_mass2);

  Double_t MTi,MTii,optimalMT2;

  MTi=calculateMT(optimal_invis_part1,s1);
  MTii=calculateMT(optimal_invis_part2,s2);
  if (MTi>=MTii){
    optimalMT2=MTi;
  }
  else if (MTi<MTii){
    optimalMT2=MTii;
  }
  else{
    optimalMT2=LONG_MAX;
  }
  return optimalMT2;
};

//this function is used to calculate alpha_t
//algorithm from Oxbridge library
void bit_increment(std::vector<Int_t>& bit_vec){
  for(Int_t bit=0; bit<( (signed)bit_vec.size() ); ++bit) {
    if(bit_vec[bit]){  //bit_vec[bit]==1
      bit_vec[bit]=0;
    }
    else{              //bit_vec[bit]==0
      bit_vec[bit]=1;
      return;
    }
  }
  // --->
  //if you got this far, the full 1 array flips over to a full 0 array
  // --->
};

//clusters jets into two pseudo-jets, with c[0] having the
//highest Et, and c[1] having the lowest Et
//YOU SHOULD CHECK TO SEE IF *ET* OR *PT* IS USED
void clusterjets(ShowerParticle* s, ShowerParticle* c, Int_t n){

  Double_t tdelET,mindelET=LONG_MAX;
  ShowerParticle s1,s2,nulls;

  if(n<=1){
    std::cout << "You can't cluster 1 jet into 2 jets, silly...\n";
    return;
  }
  else{
    std::vector<Int_t> in_clust1(n,0); //0==false 1==true
    do{
      s1=nulls;s2=nulls;
      bit_increment(in_clust1);
      for(Int_t i=0;i<n;++i){
        if(in_clust1[i]){  //we add this shower particle to cluster 1
          s1.Add(s[i]);
        }
        else{              //we add this shower particle to cluster 2
          s2.Add(s[i]);
        }
      }
      tdelET=TMath::Abs(s1.Et()-s2.Et());
      if(tdelET<mindelET){
        if(s1.Et()<s2.Et()){
          c[0]=s2;  //we arranging the cluster such that
          c[1]=s1;  // c[0] has the highest Et
        }
        else{
          c[0]=s1;
          c[1]=s2;
        }
        mindelET=tdelET;
      }
    }while (!in_clust1[n-1]);
  }
};

//Calculates the alpha_t kinematic variable
//YOU SHOULD CHECK TO SEE IF *ET* OR *PT* IS USED
//MIW MAY 2013 - Pretty sure it's ET
//definition given by arxiv 1303.2985 equation(3)
Double_t alphat(ShowerParticle* s, Int_t n){
  Double_t HT=0,HTslash,deltaHT,alphat; 
  for(Int_t i=0;i<n;++i){
    HT+=s[i].Et();
  }
  ShowerParticle smissing;
  for(Int_t i=0;i<n;++i){
    smissing.Add(s[i]);
  }
  HTslash=smissing.pt();
  ShowerParticle c[2];
  clusterjets(s,c,n);
  deltaHT=TMath::Abs( c[0].Et()-c[1].Et() );
  alphat=0.5*( (HT-deltaHT)/( TMath::Sqrt( HT*HT-HTslash*HTslash ) ) );   
//  alphat=c[1].Et()/calculateMT(c[0],c[1]); //deprecated
  return alphat;
};

#endif
