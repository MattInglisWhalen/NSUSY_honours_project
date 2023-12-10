#include <delphesv3_funcs.h>

void run_trim_delphes(char* fname){

  delphesv3* del=new delphesv3(fname);
  del->trim_delphes();

  delete del;

};
