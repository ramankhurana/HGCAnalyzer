#include "../interface/IsolationComputation.h"




std::vector<int> IsolationComputation::IsolatationRectangle(HGCalTopology* hgctopology_, 
					       std::map<HGCEEDetId,int> maxDetId,
							    int nX, int nY, int ADCCut_){
  
  std::vector<int> isodepVector;
  isodepVector.clear();

  if(maxDetId.size()==0) { std::cout<<" empty map "<<std::endl; return isodepVector;}
  else if(maxDetId.size()==1) return isodepVector;
  else {
    if(true) std::cout<<" size of map in the function = "<<maxDetId.size()<<std::endl;
    if(true) std::cout<<" det id in isolation function = "<<maxDetId.begin()->first<<std::endl;
    
    std::pair<HGCEEDetId,int> max_ = Max_Element(maxDetId);
    //maximumCell = max_.second;
    
    //if(maxDetId.size()==0)  std::cout<<" map is empty no maximum"<<std::endl;
    //else if(maxDetId.size()==1) {
    //  std::map<HGCEEDetId,int>::const_iterator iter = maxDetId.begin();
    //  maximumCell  =  iter->second;
    //}
    //else if(maxDetId.size()>1) maximumCell = max_.second;
    //else std::cout<<" map not understood "<<std::endl;
    if(true) std::cout<<" found maximum = "<<max_.first<<"   adc = "<<max_.second<<std::endl;
    for(int ixy=1;ixy<=nX;ixy++){
      std::cout<<" ixy = "<<ixy<<std::endl;
      int isodep = 0;
      for(int ix=-ixy; ix<=ixy; ix++){
	for(int iy=-ixy; iy<=ixy; iy++){
	  DetId id = hgctopology_->offsetBy(max_.first ,ix,iy);
	  HGCEEDetId newdetid = HGCEEDetId(id.rawId());
	  if(false) std::cout<<" det id of changed obj is  = "<<newdetid<<std::endl;
	  if(ix==0 && iy==0) continue ; // don't sum the seed cell
	  if(maxDetId[newdetid]>0)  std::cout<<" -------- x,y = "<<ix<<", "<<iy
					     <<" newdetid = "<<newdetid
					     <<" adc = "<<maxDetId[newdetid]<<std::endl;
	  
	  if(maxDetId.count(newdetid) >0) {
	    isodep = isodep + maxDetId[newdetid];
	    if(false) std::cout<<" passed newdetid = "<<newdetid.cell()
			       <<" adc = "<<maxDetId[newdetid]<<std::endl;
	  }
	}
      }
      isodepVector.push_back(isodep);
    }
  }
  return isodepVector;
  
}





std::vector<std::pair<int,int>> IsolationComputation::IsolatationTriangle(HGCalTopology* hgctopology_, 
									 std::map<HGCEEDetId,int> maxDetId,
									 int nX, int nY, int ADCCut_){
  
  std::vector<std::pair<int,int>> isodepVector;
  isodepVector.clear();
  if(maxDetId.size()==0) { std::cout<<" empty map "<<std::endl; return isodepVector;}
  else if(maxDetId.size()==1) return isodepVector;
  else {
    if(true) std::cout<<" size of map in the function = "<<maxDetId.size()<<std::endl;
    if(true) std::cout<<" det id in isolation function = "<<maxDetId.begin()->first<<std::endl;
    std::pair<HGCEEDetId,int> max_ = Max_Element(maxDetId);
    if(true) std::cout<<" found maximum = "<<max_.first<<"   adc = "<<max_.second<<std::endl;
    for(int ixy=1;ixy<=nX;ixy++){
      std::cout<<" ixy = "<<ixy<<std::endl;
      int isodepupper = 0;
      int isodeplower = 0;
      for(int ix=-ixy; ix<=ixy; ix++){
	for(int iy=-ixy; iy<=ixy; iy++){
	  DetId id = hgctopology_->offsetBy(max_.first ,ix,iy);
	  HGCEEDetId newdetid = HGCEEDetId(id.rawId());
	  if(false) std::cout<<" det id of changed obj is  = "<<newdetid<<std::endl;
	  if(ix==0 && iy==0) continue ; // don't sum the seed cell
	  if(maxDetId[newdetid]>0)  std::cout<<" -------- x,y = "<<ix<<", "<<iy
					     <<" newdetid = "<<newdetid
					     <<" adc = "<<maxDetId[newdetid]<<std::endl;
	  
	  if(maxDetId.count(newdetid) >0) {
	    if(ix>0 && iy>0)   isodepupper = isodepupper + maxDetId[newdetid];
	    else isodeplower = isodeplower + maxDetId[newdetid];
	    if(false) std::cout<<" passed newdetid = "<<newdetid.cell()
			       <<" adc = "<<maxDetId[newdetid]<<std::endl;
	  }
	  	  
	}
      }
      isodepVector.push_back(std::pair<int,int>(isodepupper,isodeplower));
    }
  }
  return isodepVector;
  
}




std::pair<HGCEEDetId, int> IsolationComputation::Max_Element (std::map<HGCEEDetId, int> map_){
  std::cout<<" size of map = "<<map_.size()<<std::endl;
  std::map<HGCEEDetId,int>::const_iterator iter = map_.begin();
  int maximum_= iter->second;
  HGCEEDetId   key_ = iter->first;
  //d::cout<<" maximum in function ---------"<<maximum_<<std::endl;
  if (map_.size()==0) std::cout<<" empty map can't have maximum "<<std::endl;
  else if (map_.size()==1) return std::pair<HGCEEDetId, int>(iter->first, iter->second);
  else if(map_.size()>1) {
    //d::cout<<" size of map = "<<map_.size()<<std::endl;
    for(iter=map_.begin() ; iter != map_.end(); ++iter){
      //d::cout<<" inside mal loop "<<std::endl;
      if(iter->second > maximum_){
        maximum_ = iter->second;
        key_     = iter->first;
	//d::cout<< " finding maximum = "<< key_<<"  "<<maximum_<<std::endl;
      }
    }
  }
  return std::pair<HGCEEDetId, int>(key_,maximum_);
}

