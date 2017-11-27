#include "tesselate_debug.h"


// ----------------------------------------------------------------------------  
// @@ Debug purposes
void store_points_and_curve(const ParamCurve& pc,
			    double* t,
			    unsigned int num_points,
			    const string& filename)
// ----------------------------------------------------------------------------    
{
  
  ofstream os(filename);
  pc.writeStandardHeader(os);  pc.write(os);
  vector<double> coords;
  for (double* pt = t; pt != t+num_points; ++pt) {
    Point cur_point;
    pc.point(cur_point, *pt);
    coords.push_back(cur_point[0]);
    coords.push_back(cur_point[1]);
    coords.push_back(cur_point[2]);
  }
  PointCloud<3> pcloud(coords.begin(), (int)coords.size()/3);
  pcloud.writeStandardHeader(os);
  pcloud.write(os);
  os.close();
}
