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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <float.h>
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/Body.h"

using namespace Go;
using namespace std;


int main( int argc, char* argv[] )
{
  //GoTools::init();

  if (argc < 4 )
    {
      // scale: multiplication factor used for change of units, same in all directions.
      // shift_x/y/z: translation applied after the application of the scale factor, used to realign to the origin.
      // bbmin/max_x/y/z: bounding box after scaling and translation. Only points inside the bounding box will be written out.
      cout << "Usage:  " << argv[0] << " gid_result.res gid_mesh.msh output.txt nominal_geom.g2 [scale] (optional) [shift_x shift_y shift_z] (optional) "
			            << "[bbmin_x bbmax_x bbmin_y bbmax_y bbmin_z bbmax_z] (optional)" << endl;
      return 1;
    }
  
  ifstream in_surf(argv[4]);
  CompositeModelFactory cmf(1.0e-5, 1.0e-6, 1.0e-5, 1.0e-5,  1.0e-6);
  CompositeModel* cm = cmf.createFromG2(in_surf);
  SurfaceModel* sm = dynamic_cast<SurfaceModel*>(cm); 
  shared_ptr<SurfaceModel> smm(sm);
  Body body(smm);

  vector<float> results;  
  vector<int> boundary_indices;

  ifstream in_res(argv[1]);

  std::string line;
  // Don't need first 3 lines
  std::getline(in_res, line);
  std::getline(in_res, line);
  std::getline(in_res, line);

  int index;
  float a,b,c;
  while (std::getline(in_res, line))
  { 
    if (line.substr(0,10) == string("End Values")) { break;} // eof
	std::istringstream ss(line);
	ss >> index >> a >> b >> c;
	boundary_indices.push_back(index);
	results.push_back(a);
	results.push_back(b);
	results.push_back(c);
  }

  in_res.close();

  float scale = (argc > 5) ? atof(argv[5]): 1.0f;
  float bbmin_x = (argc > 9) ? atof(argv[9]) : -FLT_MAX;
  float bbmax_x = (argc > 10) ? atof(argv[10]) : FLT_MAX;
  float bbmin_y = (argc > 11) ? atof(argv[11]) : -FLT_MAX;
  float bbmax_y = (argc > 12) ? atof(argv[12]) : FLT_MAX;
  float bbmin_z = (argc > 13) ? atof(argv[13]) : -FLT_MAX;
  float bbmax_z = (argc > 14) ? atof(argv[14]) : FLT_MAX;

  ifstream in_mesh(argv[2]);
  vector<float> pts;
  // Ignore first two lines
  std::getline(in_mesh, line);
  std::getline(in_mesh, line);
  int jndex = 0;
  int vp_count = 0;
 
 while (std::getline(in_mesh, line))
  {
    if (line.substr(0,15) == "End Coordinates") { break; } // we don't need any more
      std::istringstream ss(line);
      ss >> index;
      if (index != boundary_indices[jndex]) continue;
      ss >> a >> b >> c;
      a *= scale;
      b *= scale;
      c *= scale;

      Go::Point pt(a,b,c), clpt;
      int idx;
      double clo_par[2];
      double dist;
      sm->closestPoint(pt, clpt, idx, clo_par,dist);
      if (body.isInside(pt)) dist*=-1.0;
      //cout << "dist " << dist << endl;   
      if (dist < 0.0) { 
        vp_count++; 
        boundary_indices.erase(boundary_indices.begin()+jndex);
        results.erase(results.begin()+3*jndex,results.begin()+3*jndex+3);
        continue; 
      }  
       
      pts.push_back(a);
      pts.push_back(b);
      pts.push_back(c);
      jndex++;
  }
  in_mesh.close();

  if (pts.size() == results.size()) std::cout << "files seem to correspond ok" << std::endl;
  else std::cout << "something seems to be wrong with the input files" << std::endl;

  std::cout << "Number of points/results is " << pts.size() << " / " << results.size() << "\n\n" << std::endl; 
  std::cout << "Number of points removed (because they belong to voids) is " << vp_count << std::endl;
 
  std::vector<float> distorted;

  vector<float> shift(3);
  shift[0] = (argc > 6) ? atof(argv[6]) : 0.0f;
  shift[1] = (argc > 7) ? atof(argv[7]) : 0.0f;
  shift[2] = (argc > 8) ? atof(argv[8]) : 0.0f; 

  for (int ix=0; ix!=pts.size()/3; ++ix) {
    float a = pts[3*ix]+results[3*ix]+shift[0];
    float b = pts[3*ix+1]+results[3*ix+1]+shift[1];
    float c = pts[3*ix+2]+results[3*ix+2]+shift[2];
    if (a >= bbmin_x && a <= bbmax_x && b >= bbmin_y && b <= bbmax_y && c >= bbmin_z && c <= bbmax_z) {
      distorted.push_back(a); 
      distorted.push_back(b); 
      distorted.push_back(c); 
    }
  }

  ofstream ofs(argv[3]);
  //float barya=0.0f;
  //float baryb=0.0f;
  //float baryc=0.0f;
  int numpts = distorted.size()/3;
  for (int ix=0; ix!=numpts; ++ix) {
    //barya += pts[3*ix]/numpts;
	//baryb += pts[3*ix+1]/numpts;
	//baryc += pts[3*ix+2]/numpts;
    ofs << distorted[3*ix] << " " << distorted[3*ix+1] << " " << distorted[3*ix+2] << "\n";   
  }
  ofs << endl;
  //cout << "barycenter " << barya << " " << baryb << " " << baryc << endl; 
}
