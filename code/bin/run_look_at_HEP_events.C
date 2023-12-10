#include <delphesv3_funcs.h>

void run_look_at_HEP_events(char* fname){

  delphesv3* del=new delphesv3(fname);
  del->look_at_HEP_events();

  delete del;

};
