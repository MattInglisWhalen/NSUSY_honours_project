#ifndef MIW_ANALYSIS_H
#define MIW_ANALYSIS_H

#include "MIW_ShowerParticle.h"

template <class T>
void swap(T* ar,Int_t index1,Int_t index2){
  T temp;
  temp = ar[index1];
  ar[index1] = ar[index2];
  ar[index2] = temp;
};

//method to remove an element from an array
template <class T>
void remove(T* array, Int_t index_to_remove, Int_t& n_in_array){
  if(index_to_remove<0){
    std::cout << "MIW_analysis.h - remove() - ERROR: negative index to remove\n\n";
  }
  else{
    for(Int_t i=index_to_remove;i<n_in_array-1;++i){
      array[i]=array[i+1];
    }
    --n_in_array;
  }
};

//Simple sorting algorithm
void sort_by_pt(ShowerParticle* s,Int_t size){
  Int_t j;
  for(Int_t i=0;i<size;++i){
    j=i-1;
    if(j>=0){
      while( (s[j+1].pt()>s[j].pt()) && (j>=0) ){
        swap(s,j,j+1);
        --j;
      }
    }
  }
};

void reconstruct(ShowerParticle* s,Int_t& n_in_array,Double_t minpt=0., Double_t maxeta=LONG_MAX){
  for(Int_t i=0;i<n_in_array;++i){
    if( (s[i].pt()<minpt)||(TMath::Abs(s[i].eta())>maxeta) ){
      remove(s,i,n_in_array);
      --i;
    }
  }
};

// remove from list 1 events that are too close to list 2
void remove_overlap(ShowerParticle* s1,ShowerParticle* s2, Int_t& n1, Int_t& n2, Double_t delR,Int_t opt=0){
  Int_t removeMe;
  for(Int_t i=0;i<n1;++i){
    removeMe=0; 
    for(Int_t j=0;j<n2;++j){
      if(opt==0){  //Generic ATLAS isolation
        if(s1[i].DeltaR(s2[j])<delR){
          removeMe=1;
        }
      }
      else if(opt==1){  //ATLAS 2013 007 uses rapidity rather than eta
        if(s1[i].DeltaR_Rap(s2[j])<delR){
          removeMe=1;
        }
      }
      else if(opt==2){  //CMS 2013 008 has some weird lepton isolation with bjets
        if(s1[i].DeltaR(s2[j])<delR){
          if(s1[i].btag()>0){
            s1[i].Add(s2[j]);
            remove(s2,j,n2);
            --j;
          }
          else{
            removeMe=1;
          }
        }
      }
    }
    if(removeMe){
      remove(s1,i,n1);
      --i;
    }
  }
};

Int_t collect_leptons(ShowerParticle* e,ShowerParticle* u,ShowerParticle* l,Int_t ne,Int_t nu){
  Int_t nl=0;
  for(Int_t i=0;i<ne;++i){
    l[nl]=e[i];
    ++nl;
  }
  for(Int_t i=0;i<nu;++i){
    l[nl]=u[i];
    ++nl;
  }
  sort_by_pt(l,nl);
  return nl;
};

Int_t collect_bjets(ShowerParticle* j,ShowerParticle* b,Int_t nj,Double_t minpt=0, Double_t maxeta=LONG_MAX){
  Int_t nb=0;
  for(Int_t i=0;i<nj;++i){
    if( (j[i].btag()>0) && (j[i].pt()>minpt) && (TMath::Abs(j[i].eta())<maxeta ) ){
      b[nb]=j[i];
      ++nb;
    }
  }
  return nb;
};

Int_t collect_nonbjets(ShowerParticle* j,ShowerParticle* nonb,Int_t nj,Double_t minpt=0, Double_t maxeta=LONG_MAX){
  Int_t nnonb=0;
  for(Int_t i=0;i<nj;++i){
    if( (j[i].btag()<0.01) && (j[i].pt()>minpt) && (TMath::Abs(j[i].eta())<maxeta ) ){
      nonb[nnonb]=j[i];
      ++nnonb;
    }
  }
  return nnonb;
};


Int_t collect_taus(ShowerParticle* j,ShowerParticle* t,Int_t nj,Double_t minpt=0, Double_t maxeta=LONG_MAX){
  Int_t nt=0;
  for(Int_t i=0;i<nj;++i){
    if( (j[i].ttag()>0) && (j[i].pt()>minpt) && (TMath::Abs(j[i].eta())<maxeta ) ){
      t[nt]=j[i];
      ++nt;
    }
  }
  return nt;
};

Double_t getHT(ShowerParticle* s,Int_t njet,Int_t end=INT_MAX,Int_t start=0,Int_t opt=0){
  Double_t HT=0.;
  for(Int_t i=start;(i<end)&&(i<njet);++i){
    if(opt==0){
      HT+=s[i].pt();
    }
    else if(opt==1){
      HT+=s[i].Et();
    }
  }
  return HT;
}

//return the index of the list of shower particles sv
//correspond to the nearest particle to s1
Int_t FindNearest(ShowerParticle s1, ShowerParticle* sv, Int_t nv){
  Int_t tempnearest=-1,tempdelR=INT_MAX;
  for(Int_t i=0;i<nv;++i){
    if(s1.DeltaR(sv[i])<tempdelR){
      tempnearest=i;
    }
  }
  return tempnearest;
};

void cluster_antikt(ShowerParticle* j /*jets*/, Int_t& nj, ShowerParticle* t /*hadronic tracks*/, Int_t nt){
  // The algorithm and naming conventions are those used in http://arxiv.org/pdf/0802.1189v2.pdf

  Double_t R=0.4;     // yes, I'm hard-coding the distance parameter

  // populate the vector of pseudojets (beginning as particle tracks) that will be distributed among jet clusters
  vector<ShowerParticle> pseudo (nt,t[0]);
  for (Int_t i=0;i<nt;++i){
    pseudo[i]=t[i];
  }

  Int_t iB,ii,jj;              // indices of hadronic tracks/pseudojets to cluster together
  Double_t d_B,d_ij;           // distance measures
  Double_t tempd_B, tempd_ij;  //temporary distance measures
  Double_t fac1,fac2;          // useful variables

  //find the minimum distance measures d_B and d_ij
  while( (unsigned)pseudo.size() > 0 ){
    if( nj>300 ){
      cout << "Maximum number of jets in a single event exceeded, exiting...\n";
      exit(-1);
    } 
    d_B=d_ij=LONG_MAX;
    for(Int_t i=0;i<(unsigned)pseudo.size()-1;++i){
      tempd_B=TMath::Power( pseudo[i].pt() , -2 );
      if( tempd_B < d_B ){
        d_B=tempd_B;
        iB=i;            //found a smaller value for d_B, note the index
      }
      for(Int_t j=i+1;j<(unsigned)pseudo.size();++j){
        fac1 = TMath::Min( TMath::Power( pseudo[i].pt() , -2. ) , TMath::Power( pseudo[j].pt() , -2. ) );
        fac2 = TMath::Power( pseudo[i].DeltaR_Rap(pseudo[j])/R , 2. );
        tempd_ij=fac1*fac2;
        if ( tempd_ij < d_ij ){
          d_ij=tempd_ij;
          ii=i; jj=j;    //found a smaller value for d_ij, note the indices
        }
      }
    }
    if( d_B <= d_ij ){   // d_B == d_min, so call pseudo[iB] a jet, and remove from pseudo list
      j[nj++]=pseudo[iB];
      pseudo.erase( pseudo.begin()+iB );
      //remove(pseudo,iB,nt);
    }
    else{               // d_ij == d_min, so merge the two and remove pseudo[i] from pseudo list
      pseudo[jj]=pseudo[jj]+pseudo[ii];
      pseudo.erase( pseudo.begin()+ii );
    }
  }
};

void progress_bar(Int_t ievents,Int_t nevents){
  Int_t rat=100000*ievents/nevents;
  if(ievents<1){std::cout<<" BEGIN ["<<std::flush;}
  if(rat%1000==0){
    std::cout << '#' << std::flush;
  }
  if(ievents==(nevents-1)){std::cout<<"] DONE\n";}
};

#endif
