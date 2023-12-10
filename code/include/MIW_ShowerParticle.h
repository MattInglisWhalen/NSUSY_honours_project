#ifndef MIW_SHOWERPARTICLE_H
#define MIW_SHOWERPARTICLE_H

//C++ header files
#include <iostream> 

//root header files
#include <TMath.h>

//particle types
#define UNINITIALIZED -1
#define PHOTON 0
#define ELECTRON 1
#define MUON 2
#define TAUON 3
#define JET 4
#define MISSING 6
#define TRACK 7
#define CLUSTER 8
#define INVISIBLE 9

#define MASS_PHOTON 0.000000
#define MASS_ELECTRON 0.000511
#define MASS_MUON 0.105658
#define MASS_TAUON 1.776820
#define MASS_MISSING 0.000000

class ShowerParticle{

 protected:

  Int_t _part_type; 
  Double_t _pt,_eta,_phi,_mass; //all particles
  Double_t _charge,_hadem,_btag,_ttag; //some particles
  Double_t _etao,_phio,_vx,_vy,_vz,_vxo,_vyo,_vzo; //tracks

 public:
  //constructors and destructors
  ShowerParticle();
 ~ShowerParticle();

  //set functions
  void SetType(Int_t);
  void SetPtEtaPhiMass(Double_t,Double_t,Double_t,Double_t);
  void SetPxPyPzE(Double_t,Double_t,Double_t,Double_t);
  void SetChargeHademBtagTtag(Double_t,Double_t,Double_t,Double_t);
  void SetEtaOPhiOVxVyVzVxoVyoVzo(Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t);

  //basic get functions
  Int_t type();
  Double_t pt();
  Double_t eta();
  Double_t phi();
  Double_t mass();
  Double_t charge();
  Double_t hadem();

  //jet get functions
  Double_t btag();
  Double_t ttag();
  Double_t ntracks();

  //track get functions
  Double_t etaout();
  Double_t phiout();
  Double_t vx();
  Double_t vy();
  Double_t vz();
  Double_t vxo();
  Double_t vyo();
  Double_t vzo();

  //helper get functions
  Double_t abseta();
  Double_t mag();
  Double_t px();
  Double_t py();
  Double_t pz();
  Double_t E();
  Double_t Et();
  Double_t mt();
  Double_t rap();
  Int_t sign();

  //useful functions
  void Add(ShowerParticle);
  void NegVecSum(ShowerParticle);
  Double_t DelPhiAbs(ShowerParticle);
  Double_t DeltaR(ShowerParticle);
  Double_t DeltaR_Rap(ShowerParticle);
  Double_t InnerProduct(ShowerParticle);
  Double_t TransverseInnerProduct(ShowerParticle);
  void clear();
  void print_info();

  //operators
  ShowerParticle& operator+(ShowerParticle&);

};

ShowerParticle::ShowerParticle(){
  _part_type=UNINITIALIZED; //-1 default,0photon,1electron,2muon,3tauon,4jet
  _pt=0.;
  _eta=0.;
  _phi=0.;
  _mass=0.;

  _charge=0.;
  _hadem=0.;
  _btag=0;
  _ttag=0;

  _etao=0.;
  _phi=0.;
  _vx=0.;
  _vy=0.;
  _vz=0.;
  _vxo=0.;
  _vyo=0.;
  _vzo=0.;

};

ShowerParticle::~ShowerParticle(){
  //destructor
};

//set functions
void ShowerParticle::SetType(Int_t t){
  _part_type=t;
  if(t==PHOTON){_mass=MASS_PHOTON;}
  if(t==ELECTRON){_mass=MASS_ELECTRON;}
  if(t==MUON){_mass=MASS_MUON;}
  if(t==TAUON){_mass=MASS_TAUON;}
  if(t==MISSING){_mass=MASS_MISSING;}
}

void ShowerParticle::SetPtEtaPhiMass(Double_t p,Double_t n,Double_t f,Double_t m){
  _pt=p;
  _eta=n;
  _phi=f;
  _mass=m;
};

//set functions
void ShowerParticle::SetPxPyPzE(Double_t pxx,Double_t pyy,Double_t pzz,Double_t e){
  Double_t magp=TMath::Sqrt(pxx*pxx+pyy*pyy+pzz*pzz);
  _pt=TMath::Sqrt(pxx*pxx+pyy*pyy);
  _eta=0.5*TMath::Log( (magp+pzz)/(magp-pzz) );
  _phi=TMath::ATan2(pyy,pxx);
  if( (e*e-magp*magp) < 0. ){  //computer rounding error can suck
    _mass=-TMath::Sqrt(-(e*e-magp*magp));
  }
  else{
    _mass=TMath::Sqrt(e*e-magp*magp);
  }
};

void ShowerParticle::SetChargeHademBtagTtag(Double_t q,Double_t h,Double_t b,Double_t t){
  _charge=q;
  _hadem=h;
  _btag=b;
  _ttag=t;
};

void ShowerParticle::SetEtaOPhiOVxVyVzVxoVyoVzo(Double_t no,Double_t fo,Double_t x,Double_t y,Double_t z,Double_t xo,Double_t yo,Double_t zo){
  _etao=no;
  _phio=fo;
  _vx=x;
  _vy=y;
  _vz=z;
  _vxo=xo;
  _vyo=yo;
  _vzo=zo;
};

//get functions
Int_t ShowerParticle::type(){return _part_type;}
Double_t ShowerParticle::pt(){return _pt;}
Double_t ShowerParticle::eta(){return _eta;}
Double_t ShowerParticle::phi(){  // range of -pi < phi < pi
  Double_t f=_phi;
  while(f>TMath::Pi()){f-=2*TMath::Pi();}
  while(f<-TMath::Pi()){f+=2*TMath::Pi();}
  return f;
}
Double_t ShowerParticle::mass(){return _mass;}

//jet get functions
Double_t ShowerParticle::charge(){return _charge;}
Double_t ShowerParticle::hadem(){return _hadem;}
Double_t ShowerParticle::btag(){return _btag;}
Double_t ShowerParticle::ttag(){return _ttag;}

//track get functions
Double_t ShowerParticle::etaout(){return _etao;}
Double_t ShowerParticle::phiout(){return _phio;}
Double_t ShowerParticle::vx(){return _vx;}
Double_t ShowerParticle::vy(){return _vy;}
Double_t ShowerParticle::vz(){return _vz;}
Double_t ShowerParticle::vxo(){return _vxo;}
Double_t ShowerParticle::vyo(){return _vyo;}
Double_t ShowerParticle::vzo(){return _vzo;}

//helper get functions
Double_t ShowerParticle::abseta(){return TMath::Abs(_eta);}
Double_t ShowerParticle::mag(){return pt()*TMath::CosH(eta());}
Double_t ShowerParticle::px(){return pt()*TMath::Cos(phi());}
Double_t ShowerParticle::py(){return pt()*TMath::Sin(phi());}
Double_t ShowerParticle::pz(){return pt()*TMath::SinH(eta());}
Double_t ShowerParticle::E(){
  return TMath::Sqrt(mag()*mag()+mass()*mass());
};
Double_t ShowerParticle::Et(){
  return TMath::Sqrt( pt()*pt()+mass()*mass() );
};
Double_t ShowerParticle::mt(){
  //This is the one that makes sense to me and matches up 
  //with what gets returned by calculateMT (in MIW_kinematic_variables.h)  -- no it doesn't
  // this definition is also not boost invariant along the beam axis
  return TMath::Sqrt( mass()*mass()-pz()*pz() );

  //according to "kinematics" by the PDG, this is the correct definition
  //return Et();  //weird as hell though
}
Double_t ShowerParticle::rap(){
  return 0.5*TMath::Log( (E()+pz()) / (E()-pz()) );
}
Int_t ShowerParticle::sign(){
  if( TMath::Abs(charge())<0.0001 ){return 0;}
  return (Int_t)( charge()/TMath::Abs(charge()) );
}

//useful functions
void ShowerParticle::Add(ShowerParticle spart){
  
  _part_type=CLUSTER;
  //function variables, do a four-vector sum
  Double_t E1=E()+spart.E();
  Double_t px1=px()+spart.px();
  Double_t py1=py()+spart.py();
  Double_t pz1=pz()+spart.pz();

  Double_t mag1=TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);

  _pt=TMath::Sqrt(px1*px1+py1*py1);
  _eta=0.5*TMath::Log( (mag1+pz1)/(mag1-pz1) );
  _phi=TMath::ATan2(py1,px1);
  if( (E1*E1-mag()*mag()) < 0. ){  //computer rounding error can suck
    _mass=-TMath::Sqrt( -( E1*E1-mag()*mag() ) );
  }
  else{
    _mass=TMath::Sqrt( E1*E1-mag()*mag() );
  }
  _charge=charge()+spart.charge();
  _hadem=TMath::Sqrt(hadem()*spart.hadem()); //geometric mean, ratios don't add
  _btag=btag()+spart.btag();
  _ttag=ttag()+spart.ttag();
};

void ShowerParticle::NegVecSum(ShowerParticle spart){
  
  _part_type=CLUSTER;
  //function variables, do a four-vector subtraction
  Double_t E1=E()-spart.E();
  Double_t px1=px()-spart.px();
  Double_t py1=py()-spart.py();
  Double_t pz1=pz()-spart.pz();

  Double_t mag1=TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);

  _pt=TMath::Sqrt(px1*px1+py1*py1);
  _eta=0.5*TMath::Log( (mag1+pz1)/(mag1-pz1) );
  _phi=TMath::ATan2(py1,px1);
  if( (E1*E1-mag()*mag()) < 0. ){  //computer rounding error can suck
    _mass=-TMath::Sqrt( -( E1*E1-mag()*mag() ) );
  }
  else{
    _mass=TMath::Sqrt( E1*E1-mag()*mag() );
  }
  _charge=charge()+spart.charge();
  _hadem=TMath::Sqrt(hadem()*spart.hadem()); //geometric mean, ratios don't add
  _btag=btag()+spart.btag();
  _ttag=ttag()+spart.ttag();
};

Double_t ShowerParticle::DelPhiAbs(ShowerParticle spart){
  Double_t delphi=phi()-spart.phi();
  while(delphi>TMath::Pi()){delphi-=2*TMath::Pi();}
  while(delphi<-TMath::Pi()){delphi+=2*TMath::Pi();} 
  return TMath::Sqrt(delphi*delphi);
};

Double_t ShowerParticle::DeltaR(ShowerParticle spart){
  Double_t delphi=this->DelPhiAbs(spart);
  Double_t deleta=eta()-spart.eta();
  return TMath::Sqrt(delphi*delphi+deleta*deleta);
};

Double_t ShowerParticle::DeltaR_Rap(ShowerParticle spart){
  Double_t delphi=this->DelPhiAbs(spart);
  Double_t delrap=rap()-spart.rap();
  return TMath::Sqrt(delphi*delphi+delrap*delrap);
};
Double_t ShowerParticle::InnerProduct(ShowerParticle spart){
  return px()*spart.px()+py()*spart.py()+pz()*spart.pz();
};

Double_t ShowerParticle::TransverseInnerProduct(ShowerParticle spart){
  return px()*spart.px()+py()*spart.py();
};

void ShowerParticle::clear(){
  ShowerParticle scratch;
  (*this)=scratch;
};

void ShowerParticle::print_info(){
  std::cout << "type: " << type()
       << "  E: " << E()
       << "  x: " << px()
       << "  y: " << py()
       << "  z: " << pz() 
       << "     pt: " << pt()
       << "  eta: " << eta()
       << "  phi: " << phi()
       << "  mass: " << mass() << std::endl;
};

ShowerParticle& ShowerParticle::operator+(ShowerParticle& q){

  ShowerParticle temp;

  Double_t px1,py1,pz1,E1;
  Double_t c1,h1,b1,t1,n1;

  px1=px()+q.px();
  py1=py()+q.py();
  pz1=pz()+q.pz();
  E1=E()+q.E();

  c1=charge()+q.charge();
  h1=TMath::Sqrt(hadem()*q.hadem()); //geometric mean, ratios don't added
  b1=btag()+q.btag();
  t1=ttag()+q.ttag();

  temp.SetType(CLUSTER);
  temp.SetPxPyPzE(px1,py1,pz1,E1);
  temp.SetChargeHademBtagTtag(c1,h1,b1,t1);

  return temp;
  
};

#endif
