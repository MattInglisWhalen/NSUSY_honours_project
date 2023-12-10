#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char** argv){

  if(argc<2){
    cout << "Usage:\n   user> extract_crossx.x infile.txt\n";
    return -1;
  }

  string* infilename=new string(argv[1]);
  string* crossxfile=new string("crossx.out");
  string* crossxerrfile=new string("crossxerr.out");

  ifstream ifs(infilename->c_str());
  if(ifs.fail()){
    cout << "Could not open " << infilename->c_str() << " for reading...\n";
    return -1;
  }
  ofstream ofx(crossxfile->c_str());
  if(ofx.fail()){
    cout << "Could not open " << crossxfile->c_str() << " for reading...\n";
    return -1;
  }  
  ofstream ofe(crossxerrfile->c_str());
  if(ofe.fail()){
    cout << "Could not open " << crossxerrfile->c_str() << " for reading...\n";
    return -1;
  }  

  int mass;
  string dum1,dum2;
  double crossx,crossxerr;

  stringstream linestr;
  string line;

  getline(ifs,line);
  while( ! ifs.eof() ){
    linestr.clear();
    for(int i=0;i<line.length();++i){
      if(line[i]=='?'){line[i]=' ';}
    }
    linestr.str(line);
    linestr>>mass>>dum1>>crossx>>dum1>>crossxerr;
    ofx<<crossx<<",";
    ofe<<crossx*crossxerr/100.<<",";
    getline(ifs,line);   
  }
  ifs.close();
  ofx.close();
  ofe.close();


  return 0;
};
