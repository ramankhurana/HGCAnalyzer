#ifndef IsolationComputation_H
#define IsolationComputation_H

#include <memory>
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include "TMath.h"
#include <vector>
#include "TH1F.h"
#include "TRandom.h"
#include "TFile.h"
#include "TString.h"
#include <sstream>
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h" 

class IsolationComputation {
 public :
  IsolationComputation(){};
  ~IsolationComputation(){};
  
  // Return the sum of ADC counts in a rectange around the maximum ADC cell
  // it will look for nX (nY) cells around the hot cell in both the directions. 
  // the adc to be considered should have atleast ADCCut ADC counts
  std::vector<int> IsolatationRectangle(HGCalTopology* hgctopology_, 
					std::map<HGCEEDetId,int> maxDetId,
					int nX, int nY, int ADCCut_);
  
  std::vector<std::pair<int,int>> IsolatationTriangle(HGCalTopology* hgctopology_, 
						      std::map<HGCEEDetId,int> maxDetId,
						      int nX, int nY, int ADCCut_);
  
  std::pair<HGCEEDetId, int> Max_Element (std::map<HGCEEDetId, int> map_);
  // class variables
  
  // class functions
  
  
 private :
  
  // private members
  
};
#endif


