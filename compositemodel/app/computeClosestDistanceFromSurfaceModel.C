/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */



#include <vector>
#include <fstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"


using namespace Go;
using namespace std;



static bool abs_comp(float a, float b)
{
    return (std::fabs(a) < std::fabs(b));
}


int main( int argc, char* argv[] )
{
  GoTools::init();

  if (argc < 4 || argc > 5)
    {
      cout << "Usage:  " << argv[0] << " surfacefile.g2 pointfile.txt outfile.txt [outdistances.txt] (optional)" << endl;
      return 1;
    }

  ifstream in_surf(argv[1]);
  CompositeModelFactory cmf(1.0e-5, 1.0e-6, 1.0e-5, 1.0e-5,  1.0e-6);
  CompositeModel* sm = cmf.createFromG2(in_surf);

  ifstream in_pts(argv[2]);
  vector<float> pts;

  while (!in_pts.eof())
    {
      for (int j = 0; j < 3; ++j)
	{
	  float f;
	  in_pts >> f;
	  pts.push_back(f);
	}
      Utils::eatwhite(in_pts);
    }
  in_pts.close();

  vector<float> closest_pts;
  for (int ix=0; ix!=pts.size()/3;++ix) {
   Go::Point pt(pts[3*ix],pts[3*ix+1],pts[3*ix+2]), clpt;
   int idx;
   double clo_par[2];
   double dist;
   sm->closestPoint(pt, clpt, idx, clo_par,dist);
   closest_pts.push_back((float) dist);
  }

  auto abs_max_it = max_element(std::begin(closest_pts), std::end(closest_pts), abs_comp);
  auto abs_min_it = min_element(std::begin(closest_pts), std::end(closest_pts), abs_comp);
  auto max_it = max_element(std::begin(closest_pts), std::end(closest_pts));
  auto min_it = min_element(std::begin(closest_pts), std::end(closest_pts));
  
  ofstream ofs(argv[3]);
  ofs << *min_it << "\n";
  ofs << *max_it << "\n\n";  
 
  ofs << pts[3*(min_it-closest_pts.begin())] << " " << pts[3*(min_it-closest_pts.begin())+1] << " " << pts[3*(min_it-closest_pts.begin())+2] << "\n";
  ofs << pts[3*(max_it-closest_pts.begin())] << " " << pts[3*(max_it-closest_pts.begin())+1] << " " << pts[3*(max_it-closest_pts.begin())+2] << "\n" << endl;

  if (argc > 4) {
	ofstream ofsdist(argv[4]);
    for (int i = 0; i < closest_pts.size(); i += 1)
      ofsdist << closest_pts[i] << "\n";
  }
   

}
