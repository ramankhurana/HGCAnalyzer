#include "../interface/hitsinfo.h"
//ClassImp(hitsinfo);


hitsinfo::hitsinfo() {
  x=y=z=time=energy=phi=eta=0.0;
  cell=sector=layer=0;
  subdetname="";
}

hitsinfo::~hitsinfo(){}
