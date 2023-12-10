#ifndef GRIDPOINT_H
#define GRIDPOINT_H

#include <iostream>

#define SBOT1  1000005
#define STOP1  1000006
#define GLUINO 1000021
#define NEUT1  1000022
#define NEUT2  1000023
#define CHRG1  1000024
#define NEUT3  1000025
#define NEUT4  1000035
#define CHRG2  1000037
#define SBOT2  2000005
#define STOP2  2000006

class Gridpoint{

 private:
  Int_t _nSR;
  Double_t*     _eff[3];  // eff[mode][SR]
  Int_t     _memflag[3];  //tells me if i've assigned memory to _eff[i]
  Double_t* _xyzvals[3];  //pointers to the masses i'm using for x, y or z
  Int_t     _xyzpart[3];  //tells me which particle mass is x, y or z
  Double_t  _mb1,_mt1,_mgo,_mn1,_mn2,_mx1,_mn3,_mn4,_mx2,_mb2,_mt2;  //this order is important

 public:
  //constructors and destructors
  Gridpoint(Int_t=1);
  Gridpoint(const Gridpoint&);
 ~Gridpoint();

  //get functions
  Int_t nSR() const;
  Double_t mpart(Int_t) const;
  Double_t x() const;
  Double_t y() const;
  Double_t z() const;
  Double_t xyzvals(Int_t) const;
  Double_t xpart() const;
  Double_t ypart() const;
  Double_t zpart() const;
  Double_t xyzpart(Int_t) const;
  Double_t eff(Int_t,Int_t) const;  //mode, SR

  //set functions
  void setx(Double_t);
  void sety(Double_t);
  void setz(Double_t);
  void setvali(Double_t,Int_t);
  void setxyzvals(Double_t,Double_t,Double_t);
  void setxyzpart(Int_t,Int_t,Int_t);
  void setnxgbt(Double_t,Double_t,Double_t,Double_t,Double_t);
  void setEff(Double_t,Int_t,Int_t);  //val, mode, SR

  //useful functions
  Double_t eff_eff(Double_t,Int_t,Int_t=0) const; //effective efficiency
  
  //operators
  const Gridpoint& operator=(const Gridpoint&);

};

Gridpoint::Gridpoint(Int_t nSR) : _nSR(nSR){
  _mb1=_mt1=_mgo=_mn1=_mn2=_mx1=_mn3=_mn4=_mx2=_mb2=_mt2=LONG_MIN;
  memset(_eff,NULL,3*sizeof(Double_t*));
  memset(_memflag,0,3*sizeof(Int_t));
  memset(_xyzvals,NULL,3*sizeof(Double_t*));
  memset(_xyzpart,-1,3*sizeof(Int_t*));
  for(Int_t i=0;i<3;++i){
    //making sure the efficiencies are initialized to zero
    _eff[i]=new Double_t[nSR];
    memset(_eff[i],0.,nSR*sizeof(Double_t));
    //letting the Gridpoint know we've created memory for _eff[i]
    _memflag[i]=1;
  }
};

Gridpoint::Gridpoint(const Gridpoint& GP) : 
  _nSR(GP._nSR),
  _mb1(GP._mb1),
  _mt1(GP._mt1),
  _mgo(GP._mgo),
  _mn1(GP._mn1),
  _mn2(GP._mn2),
  _mx1(GP._mx1),
  _mn3(GP._mn3),
  _mn4(GP._mn4),
  _mx2(GP._mx2),
  _mb2(GP._mb2),
  _mt2(GP._mt2){
  setxyzpart(GP.xpart(),GP.ypart(),GP.zpart());
  memset(_memflag,0,3*sizeof(Int_t));
  memset(_eff,NULL,3*sizeof(Double_t*));
  for(Int_t i=0;i<3;++i){
    //making memory for the efficiencies
    _eff[i]=new Double_t[_nSR];
    for(Int_t j=0;j<_nSR;++j){
      _eff[i][j]=GP.eff(i,j);
    }
    //letting the Gridpoint know we've created memory for the _eff[i] 
    _memflag[i]=1;
  }
};

Gridpoint::~Gridpoint(){
  for(Int_t i=0;i<3;++i){
    if(_memflag[i]==1){
      delete[] _eff[i];
      _memflag[i]=0;
    }
  }
};

Int_t Gridpoint::nSR() const {return _nSR;}
Double_t Gridpoint::mpart(Int_t part) const{
       if(part== SBOT1){return _mb1;}
  else if(part== STOP1){return _mt1;}
  else if(part==GLUINO){return _mgo;}
  else if(part== NEUT1){return _mn1;}
  else if(part== NEUT2){return _mn2;}
  else if(part== CHRG1){return _mx1;}
  else if(part== NEUT3){return _mn3;}
  else if(part== NEUT4){return _mn4;}
  else if(part== CHRG2){return _mx2;}
  else if(part== SBOT2){return _mb2;}
  else if(part== STOP2){return _mt2;}
  return -1;
};

Double_t Gridpoint::x() const {return this->xyzvals(0);}
Double_t Gridpoint::y() const {return this->xyzvals(1);}
Double_t Gridpoint::z() const {return this->xyzvals(2);}
Double_t Gridpoint::xyzvals(Int_t i) const {
  if(_xyzpart[i]>0){return *(_xyzvals[i]);}
  return LONG_MIN;
}
Double_t Gridpoint::xpart() const {return _xyzpart[0];}
Double_t Gridpoint::ypart() const {return _xyzpart[1];}
Double_t Gridpoint::zpart() const {return _xyzpart[2];}
Double_t Gridpoint::xyzpart(Int_t i) const {return _xyzpart[i];}
Double_t Gridpoint::eff(Int_t imode,Int_t iSR) const {return _eff[imode][iSR];}

//set functions
void Gridpoint::setx(Double_t xx){
  setvali(xx,0);
};
void Gridpoint::sety(Double_t yy){
  setvali(yy,1);
};
void Gridpoint::setz(Double_t zz){
  setvali(zz,2);
};
void Gridpoint::setvali(Double_t val,Int_t i){
  if(_xyzpart[i]>0){*(_xyzvals[i])=val;}
};
void Gridpoint::setxyzvals(Double_t mx,Double_t my,Double_t mz){
  this->setx(mx);
  this->sety(my);
  this->setz(mz);
};
//set functions
void Gridpoint::setxyzpart(Int_t whichx,Int_t whichy,Int_t whichz){
  _xyzpart[0]=whichx;
  if(whichx==whichy){return;}
  _xyzpart[1]=whichy;
  if( (whichx==whichz)||(whichy==whichz) ){return;}
  _xyzpart[2]=whichz;
  for(Int_t ixyz=0;ixyz<3;++ixyz){
         if(_xyzpart[ixyz]==SBOT1) {_xyzvals[ixyz]=&_mb1;}
    else if(_xyzpart[ixyz]==STOP1) {_xyzvals[ixyz]=&_mt1;}
    else if(_xyzpart[ixyz]==GLUINO){_xyzvals[ixyz]=&_mgo;}
    else if(_xyzpart[ixyz]==NEUT1) {_xyzvals[ixyz]=&_mn1;}
    else if(_xyzpart[ixyz]==NEUT2) {_xyzvals[ixyz]=&_mn2;}
    else if(_xyzpart[ixyz]==CHRG1) {_xyzvals[ixyz]=&_mx1;}
    else if(_xyzpart[ixyz]==NEUT3) {_xyzvals[ixyz]=&_mn3;}
    else if(_xyzpart[ixyz]==NEUT4) {_xyzvals[ixyz]=&_mn4;}
    else if(_xyzpart[ixyz]==CHRG2) {_xyzvals[ixyz]=&_mx2;}
    else if(_xyzpart[ixyz]==SBOT2) {_xyzvals[ixyz]=&_mb2;}
    else if(_xyzpart[ixyz]==STOP2) {_xyzvals[ixyz]=&_mt2;}
  }
};
void Gridpoint::setnxgbt(Double_t mn,Double_t mx,Double_t mg,Double_t mb,Double_t mt){
  _mn1=mn;
  _mx1=mx;
  _mgo=mg;
  _mb1=mb;
  _mt1=mt;
};
void Gridpoint::setEff(Double_t val,Int_t imode,Int_t iSR){_eff[imode][iSR]=val;}

const Gridpoint& Gridpoint::operator=(const Gridpoint& GP){
  _nSR=GP.nSR();
  _mb1=GP.mpart(SBOT1);
  _mt1=GP.mpart(STOP1);
  _mgo=GP.mpart(GLUINO);
  _mn1=GP.mpart(NEUT1);
  _mn2=GP.mpart(NEUT2);
  _mx1=GP.mpart(CHRG1);
  _mn3=GP.mpart(NEUT3);
  _mn4=GP.mpart(NEUT4);
  _mx2=GP.mpart(CHRG2);
  _mb2=GP.mpart(SBOT2);
  _mt2=GP.mpart(STOP2);
  setxyzpart(GP.xpart(),GP.ypart(),GP.zpart());
  for(Int_t i=0;i<3;++i){
    if(_memflag[i]==1){
      delete[] _eff[i];
      _memflag[i]=0;
    }
    _eff[i]=new Double_t[_nSR];
    _memflag[i]=1;
    for(Int_t j=0;j<_nSR;++j){
      _eff[i][j]=GP.eff(i,j);
    }
  }
  return GP;
};

//BRT refers to the probability that the particle will go through the top decay mode
Double_t Gridpoint::eff_eff(Double_t BRT, Int_t SR,Int_t opt) const {
  Double_t effeff;

  if(opt==0){   // the real definition if we have all 3 mode efficiencies
    effeff = _eff[0][SR]*BRT*BRT + 2*_eff[1][SR]*BRT*(1-BRT) + _eff[2][SR]*(1-BRT)*(1-BRT);
  } else if(opt==1){  // "Rescaling" dominant t-topology
    effeff = _eff[0][SR]*BRT*BRT;  
  } else if(opt==2){  // "Rescaling" dominant b-topology
    effeff = _eff[2][SR]*(1-BRT)*(1-BRT);  
  } else if(opt==3){  // "Rescaling" both topologies
    effeff = _eff[0][SR]*BRT*BRT + _eff[2][SR]*(1-BRT)*(1-BRT);
  } else if(opt==4){  //educated guess for what it should be... cross mode efficiency is a
                      // geometric mean between tmode and bmode efficiencies
    effeff = _eff[0][SR]*BRT*BRT + 2*TMath::Sqrt(_eff[0][SR]*_eff[2][SR])*BRT*(1-BRT) + _eff[2][SR]*(1-BRT)*(1-BRT);  
  }
  else{
    effeff=0.;
  }
  return effeff;
};

#endif