#include <delphesv3_funcs.h>

void run_root_analysis(char* ansystag){

  delphesv3* del=new delphesv3(ansystag);
  del->RunAnalysis();

  std::cout << endl;  //for formatting

  delete del;

};

