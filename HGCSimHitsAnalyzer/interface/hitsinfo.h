#ifndef hitsinfo_h
#define hitsinfo_h
#include <iostream>
#include "TObject.h" 
using namespace std;
class hitsinfo : public TObject{
 public :
  hitsinfo() ;
  ~hitsinfo() ;
  
  double x, y, z, time, energy, phi, eta ;
  int    cell, sector, layer;
  std::string subdetname;
 private :
  ClassDef(hitsinfo,1);
  
};

#endif
