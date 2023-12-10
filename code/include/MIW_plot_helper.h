#ifndef PLOT_HELPER_HH
#define PLOT_HELPER_HH

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include <MIW_useful_functions.h>
#include <MIW_search_data.h>
#include <MIW_Gridpoint.h>

#define CONF_LEVEL 0.95
#define TOPMASS 172.9
#define BOTMASS 4.44

using namespace std;

typedef vector< Gridpoint > vec1;
typedef vector<    vec1   > vec2;
typedef vector<    vec2   > vec3;

//global variable holding all the gridpoint data
Int_t nSR;
Double_t *data,*smpred,*error_smpred,lumin,error_lumin,*plot_opts,*UL95;
vec1 GP;  //base gridpoints
vec1 exclGP,exclGP_tmode,exclGP_bmode,exclGP_both,exclGP_geom;  //excluded gridpoints
vec1 ordGP,ordGP_tmode,ordGP_bmode,ordGP_both,ordGP_geom;  //ordered gridpoints to make a nice contour

void get_search_data(char* ansystag){
  nSR=get_nSR(ansystag);
  data=get_data(ansystag);
  smpred=get_smpred(ansystag);
  error_smpred=get_error_smpred(ansystag);
  lumin=get_lumin(ansystag);
  error_lumin=get_error_lumin(ansystag);
  plot_opts=get_plot_opts(ansystag);
  UL95=get_UL95(ansystag);
};

//I don't need masses to be the same down to 5 decimals, I only need integer accuracy
//Scratch that, it depends on how large the values are. It should be done in relative error
//ie tolerance is in percent error (0-100)
Int_t about_the_same_rel(Double_t x, Double_t y, Double_t tol){
  Double_t numerator=40000*(x-y)*(x-y);
  Double_t denominator=(x+y)*(x+y);
  if(numerator<0.001){return 1;}
  if( numerator/denominator<(tol)*(tol) ){return 1;}
  return 0;
};

Int_t about_the_same_abs(Double_t x,Double_t y,Double_t diff){
  if ( (x-y)*(x-y)<diff*diff ){ return 1;}
  return 0;
};

//so we don't create two gridpoints for different modes at the same location
Int_t findGP(Gridpoint gp){
  for(Int_t i=0;i<(signed)GP.size();++i){
    if(  about_the_same_rel(GP[i].x(),gp.x(),1) 
      && about_the_same_rel(GP[i].y(),gp.y(),1) 
      && about_the_same_rel(GP[i].z(),gp.z(),1)  ){
      return i;
    }
  }
  return -1;
};

//read ansystag.eff held in a directory dictated by mode
Int_t read_data(char* stortag, char* ansystag, char* ztag, Int_t mode,Int_t xopt, Int_t yopt, Int_t zopt){ //mode: 0==TT, 1==TB, 2==BB

  Double_t tmu,tM1,tM2,tM3,tMQ3,tMU3,tN1mass,tX1mass,tGmass,tB1mass,tT1mass,teff;  //temporary variables
  Gridpoint tGP(nSR);
  tGP.setxyzpart(xopt,yopt,zopt);
  Int_t gpind;

  //go to the analysis director to do the reading in
  string* effname=new string(stortag);
  if(mode==0){effname->append("_TMODE/");}
  else if(mode==1){effname->append("_CMODE/");}
  else if(mode==2){effname->append("_BMODE/");}
  else{cout << "read_data: Invalid mode (" << mode << "), try again.\n";exit(-1);}
  effname->append(ansystag); effname->append("/");
  if ( yopt != NEUT1 ) {effname->append(ztag);effname->append("_");}
  effname->append(ansystag);  effname->append(".eff");

  ifstream ff;
  string* line=new string("scratch"); line->clear();
  stringstream linestr;

  //get efficiencies
  ff.open(effname->c_str(),ios::in);
  if( !ff.is_open() ){
    cout << "Could not open " << *effname << " for reading.\n";
    exit(-1);
  }
  getline(ff,(*line));  // tag
  getline(ff,(*line));  // format
  //assuming the format of mu, M1, M2, M3, MQ3, MU3, N1 Mass, X1 Mass, G Mass, B1Mass, T1 Mass
  while ( !ff.eof() ){
    line->clear();
    getline(ff,(*line));
    for(Int_t i=0;i<(signed)line->length();++i){
      if((*line)[i]==','){
        (*line)[i]=' ';
      }
    }
    linestr.clear();
    linestr.str((*line));
    linestr>>tmu>>tM1>>tM2>>tM3>>tMQ3>>tMU3>>tN1mass>>tX1mass>>tGmass>>tB1mass>>tT1mass;
    tGP.setnxgbt(tN1mass,tX1mass,tGmass,tB1mass,tT1mass);
    gpind=findGP(tGP);
    for(Int_t iSR=0;iSR<nSR;++iSR){
      linestr>>teff;
      if(gpind<0){tGP.setEff(teff,mode,iSR);}
      else{GP[gpind].setEff(teff,mode,iSR);}
    }
    if(gpind<0){GP.push_back(tGP);}
  }
  ff.close();

  delete effname;
  delete line;

  return 0;  
};

void swapGP(vector<Gridpoint>& gp, Int_t index1,Int_t index2){
  Gridpoint temp;
  temp = gp[index1];
  gp[index1] = gp[index2];
  gp[index2] = temp;
};

//constant x, then y, then z
void sort_xyz(){
  Int_t j;
  for(Int_t i=1;i<(signed)GP.size();++i){
    j=i-1;
    while(j>=0){
      if( about_the_same_abs(GP[j+1].x(),GP[j].x(),5) ){
        if( about_the_same_abs(GP[j+1].y(),GP[j].y(),5) ){
          if ( GP[j+1].z()<GP[j].z() ){swapGP(GP,j,j+1);}
          else { /* do nothing */ ; }
        }
        else if( GP[j+1].y()<GP[j].y() ){swapGP(GP,j,j+1);}
        else { /* do nothing */ ; }
      }
      else if( GP[j+1].x()<GP[j].x() ){swapGP(GP,j,j+1);}
      else{ /* do nothing */ ; }
      --j;
    }
  }
};

//constant z, then y, then x
void sort_zyx(){
  Int_t j;
  for(Int_t i=1;i<(signed)GP.size();++i){
    j=i-1;
    while(j>=0){
      if( about_the_same_rel(GP[j+1].z(),GP[j].z(),1) ){
        if( about_the_same_rel(GP[j+1].y(),GP[j].y(),1) ){
          if ( GP[j+1].x()<GP[j].x() ){swapGP(GP,j,j+1);}
          else { /* do nothing */ ; }
        }
        else if( GP[j+1].y()<GP[j].y() ){swapGP(GP,j,j+1);}
        else { /* do nothing */ ; }
      }
      else if( GP[j+1].z()<GP[j].z() ){swapGP(GP,j,j+1);}
      else{ /* do nothing */ ; }
      --j;
    }
  }
};

void permute_GP_axes(){
  Int_t tempxpart,tempypart,tempzpart;
  for(Int_t i=0;i<(signed)GP.size();++i){
    tempxpart=GP[i].xpart();
    tempypart=GP[i].ypart();
    tempzpart=GP[i].zpart();
    GP[i].setxyzpart(tempzpart,tempxpart,tempypart);
  }
};

//I have to assume that the gridpoint vector has been sorted
//refine then goes along the z-axis (constant x and y) 
// and inserts gridpoints with interpolated efficiencies
void refine(){
  
  Gridpoint newGP;

  for(Int_t ipermute=0;ipermute<3;++ipermute){
    sort_xyz();
    for(Int_t i=1;i<(signed)GP.size();++i){
      if( about_the_same_rel(GP[i].x(),GP[i-1].x(),0.5)
       && about_the_same_rel(GP[i].y(),GP[i-1].y(),0.5) ) {
        newGP=GP[i];
        newGP.setxyzvals( (GP[i].x()+GP[i-1].x())/2. ,
                          (GP[i].y()+GP[i-1].y())/2. ,
                          (GP[i].z()+GP[i-1].z())/2. );
        for(Int_t mode=0;mode<3;++mode){
          for(Int_t iSR=0;iSR<newGP.nSR();++iSR){
            newGP.setEff( (GP[i].eff(mode,iSR)+GP[i-1].eff(mode,iSR))/2. , mode, iSR );
          }
        }
        GP.insert(GP.begin()+i,newGP);  //the only thing that changes
        ++i;
      }
    }
    permute_GP_axes();
  }
  sort_xyz();
};

//naive approach to CLs - exact for only one signal region
Double_t CL_s(Double_t BR, Int_t SR, Int_t iGP, Int_t proc,Int_t opt=1){  //proc=1 for stop, proc=2 for gluino pair production
                                                                          //opt=0 for exact, 1 for tmode, 2 for b mode, 3 for both, 4 for geo,etric
  Double_t effeff=GP[iGP].eff_eff(BR,SR,opt);
  Double_t nevsignal=0.;
  if(proc==1){nevsignal=crossx(GP[iGP].mpart(STOP1),proc)*lumin*effeff;}
  if(proc==2){nevsignal=crossx(GP[iGP].mpart(GLUINO),proc)*lumin*effeff;}
  Double_t nuSUSY =nevsignal+smpred[SR];
  Double_t nuSM  = smpred[SR];

  return ( ROOT::Math::poisson_cdf( (Int_t)data[SR],nuSUSY )/ROOT::Math::poisson_cdf( (Int_t)data[SR],nuSM ) );  
};

Double_t p_value(Double_t BR, Int_t iGP,Int_t proc,Int_t opt=0){

  //p-value for mass spectrum point is the max p-value across all the SR's
  // calculated using the CLs prescription
  Double_t max_pval=0.,temp_pval;
  for(Int_t iSR=0;iSR<nSR;++iSR){
    temp_pval=1-CL_s(BR,iSR,iGP,proc,opt);
    if(temp_pval>max_pval){
      max_pval=temp_pval;
    }
  }
  return max_pval;
};

Int_t sigAe_greater_than_UL(Double_t BR, Int_t iGP,Int_t proc,Int_t opt=0){

  Double_t effeff,sigAe; 
  for(Int_t iSR=0;iSR<nSR;++iSR){
    effeff=GP[iGP].eff_eff(BR,iSR,opt);
    sigAe=crossx(GP[iGP].mpart(GLUINO),proc)*effeff*1000.;  //1000 is to go from pb to fb
    if(sigAe>=UL95[iSR]){
      return 1;
    }
  }
  return 0;
};

//this function is specific to the go > ttxn1/bbxn1 decay process
Int_t allowed_go_ttxbbxn1(Gridpoint grid, Double_t br){
  if( (br>0.) && (grid.mpart(GLUINO)<(2*TOPMASS+grid.mpart(NEUT1))) ){
    return 0;
  }
  if( (br<1.) && (grid.mpart(GLUINO)<(2*BOTMASS+grid.mpart(NEUT1))) ){
    return 0;
  }
  return 1;
};
//this function is specific to the go > t1 t~ decay process
Int_t allowed_go_t1tx(Gridpoint grid){
  if( grid.mpart(GLUINO)<(TOPMASS+grid.mpart(STOP1)) ){
    return 0;
  }
  return 1;
};
//this function is specific to the t1 > t n1 decay process
Int_t allowed_t1_tn1_bx1(Gridpoint grid, Double_t br){
  if( (br>0.) && (grid.mpart(STOP1)<(TOPMASS+grid.mpart(NEUT1))) ){
    return 0;
  }
  if( (br<1.) && (grid.mpart(STOP1)<(BOTMASS+grid.mpart(CHRG1))) ){
    return 0;
  }
  return 1;
};
Int_t allowed(Gridpoint grid,Double_t br,Int_t dec1,Int_t dec2){
  if(dec1==1){return allowed_go_ttxbbxn1(grid,br);}
  if( (dec1==2)&&(dec2==3)&&(allowed_go_t1tx(grid)) ){return allowed_t1_tn1_bx1(grid,br);}
  return 0;
};
void find_excluded_GPs(Double_t BR, Int_t proc,Int_t DEC1,Int_t DEC2=0,Int_t opt=0){
  for(Int_t i=0;i<(signed)GP.size();++i){
    //calculate the p-value (using the CL_s 
    // statistic) for each gridpoint 
    if(p_value(BR,i,proc,opt)>=CONF_LEVEL){
      if(allowed(GP[i],BR,DEC1,DEC2)){
        if(opt==0){exclGP.push_back(GP[i]);}
        else if(opt==1){exclGP_tmode.push_back(GP[i]);}
        else if(opt==2){exclGP_bmode.push_back(GP[i]);}
        else if(opt==3){exclGP_both.push_back(GP[i]);}
        else if(opt==4){exclGP_geom.push_back(GP[i]);}
      }
    }
  }
};

void find_excluded_GPs_UL95(Double_t BR, Int_t proc,Int_t DEC1,Int_t DEC2=0,Int_t opt=0){
  for(Int_t i=0;i<(signed)GP.size();++i){
    if(sigAe_greater_than_UL(BR,i,proc,opt)){
      if(allowed(GP[i],BR,DEC1,DEC2)){
        if(opt==0){exclGP.push_back(GP[i]);}
        if(opt==1){exclGP_tmode.push_back(GP[i]);}
        if(opt==2){exclGP_bmode.push_back(GP[i]);}
        if(opt==3){exclGP_both.push_back(GP[i]);}
        if(opt==4){exclGP_geom.push_back(GP[i]);}      }
    }
  }
}

//for converting a mass to a string
string mass_string(Double_t mass, char* front=NULL, char* back=NULL){
  string ret_string;
  if(front){
    ret_string=front;
  }
  mass=TMath::Ceil(mass);
  if(mass>=1000){
    ret_string+=( ((Int_t)TMath::Floor(mass/1000))%10 + 48);
  }
  if(mass>=100){
    ret_string+=( ((Int_t)TMath::Floor(mass/100))%10 + 48);
  }
  if(mass>=10){
    ret_string+=( ((Int_t)TMath::Floor(mass/10))%10 + 48);
  }
  ret_string+=( ((Int_t)TMath::Floor(mass/1))%10 + 48 );
  if(back){
    ret_string+=back;
  }
  return ret_string;
};

//for converting a branching ratio to a string
// returns to a precision of 0.01
string BR_string(Double_t BR, char* front=NULL, char* back=NULL){
  string ret_string;
  if(front){
    ret_string=front;
  }
  if(BR>=1.){ret_string+=( ((Int_t)TMath::Floor(BR))%10 + 48);ret_string+=".";}
  else{ret_string+="0.";}
  if(BR>=0.1){ret_string+=( ((Int_t)TMath::Floor(BR*10))%10 + 48);}
  else{ret_string+="0";}
  if(BR>=0.01){ret_string+=( ((Int_t)TMath::Floor(BR*100))%10 + 48);}
  else{ret_string+="0";}
  if(back){ret_string+=back;}
  return ret_string;
};

vec1 order_excluded_GPs(vec1 thisexclGP){
  vec1 ordGP;
  if(thisexclGP.size()==0){return ordGP;}
  vec1* left=new vec1;
  vec1* right=new vec1;
  left->push_back(thisexclGP[0]);
  for(Int_t i=1;i<(signed)thisexclGP.size();++i){
    if(!about_the_same_rel(thisexclGP[i-1].y(),thisexclGP[i].y(),1)){
      right->push_back(thisexclGP[i-1]);
      left->push_back(thisexclGP[i]);
    }  
  }
  right->push_back( thisexclGP[thisexclGP.size()-1] );

  //go up along the left
  for(Int_t i=0;i<(signed)left->size();++i){
    ordGP.push_back(left->at(i));
  }
  //and down along the right
  for(Int_t i=((signed)right->size())-1;i>=0;--i){
    ordGP.push_back(right->at(i));

  }
  //and close the contour
  ordGP.push_back(left->at(0));
  return ordGP;
};

//for plotting purposes
//xx[0]==Branching Ratio
//pars[0]==Mass Point #
Double_t p_value_plot(Double_t* xx, Double_t* pars){
  return p_value(xx[0],(Int_t)pars[0],(Int_t)pars[1]);
};

#endif
