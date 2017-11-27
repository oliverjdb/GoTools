#ifndef _CLIPPED_DOMAIN_2D
#define _CLIPPED_DOMAIN_2D

#include <vector>

namespace Go {

class ClippedDomain2D
{
public:
  ClippedDomain2D(const std::vector<std::vector<double>>& bounding_loops,
		  unsigned int num_boxes_x, unsigned int num_boxes_y) 
  
private:
  std::vector<std::vector<double>> loops_;
  std::vector<int> boxes_;
  
}
  
};

#endif
