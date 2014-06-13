#ifndef energysum_h
#define energysum_h
#include "TObject.h"

class energysum : public TObject{
 public :
  energysum() ;
  ~energysum();
  double e15, e25, e50, e100, e250, e1000;
 private :
  ClassDef(energysum,1);
};

#endif
