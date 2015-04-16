#include "interface/DetectorEffects.h"
#include <ctime>
#include <iostream>
using namespace std;
int main(int argc, char *argv[]){
  time_t start, end;
  time(&start);
  
  if(argc<3) {
    std::cout<< " insufficient input provided "<<std::endl;
    std::cout<< " usage : ./main.exe input.root output.root "<<std::endl;
    return 0;
  }
  
  TString input  = argv[1];
  TString output = argv[2];
  TTree *tree=0;
  DetectorEffects* detEffects;
  detEffects = new DetectorEffects(tree,input,output);
  detEffects->Loop();
  detEffects->DetectorEffects::~DetectorEffects();
  
  
  time(&end);
  std::cout<<" time used is = "<<-(start-end)<<" seconds"<<std::endl;
  
  return 0;
}
