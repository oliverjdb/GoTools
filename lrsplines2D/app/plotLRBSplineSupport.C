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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineEvalGrid.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

using namespace Go;
using namespace std;

typedef LRSplineSurface::ElementMap::const_iterator elem_it;
typedef std::vector<double>::const_iterator vecit;
typedef std::vector<LRBSpline2D*>::iterator lrb_it;

vector<double> make_regular_kvec(int degree, int num_intervals, double parmin, double parmax)// bool full_mult_bd=false)
//------------------------------------------------------------------------------
{
	assert(parmin < parmax);

	vector<double> result;
	const double interval_size = (parmax - parmin) / num_intervals;

	//if (full_mult_bd) for (int i = 0; i != degree; ++i) result.push_back(parmin);
	for (int i = 0; i != num_intervals+1; ++i) result.push_back(parmin + i * interval_size);
	//if (full_mult_bd) for (int i = 0; i != degree; ++i) result.push_back(parmax);

	return result;
}

int main(int argc, char *argv[])
{
  if (argc != 2) {
    cout << "Usage: ./LRandMS <output prefix (many files will be generated)>" << endl;
    return -1;
  }

  int deg_x = 2;
  int deg_y = 2;
  int elements_x = 5;
  int elements_y = 5;
  double parmin_x = 0;
  double parmin_y = 0;
  double parmax_x = 1.0;
  double parmax_y = 1.0;

  const vector<double> kvec_x = make_regular_kvec(deg_x, elements_x, parmin_x, parmax_x);
  const vector<double> kvec_y = make_regular_kvec(deg_y, elements_y, parmin_y, parmax_y);

  for (int ix = 0; ix != kvec_x.size(); ++ix) {
	  cout << kvec_x[ix] << " ";
  } cout << endl;

  for (int ix = 0; ix != kvec_y.size(); ++ix) {
	  cout << kvec_y[ix] << " ";
  } cout << endl;

  LRSplineSurface lrs(deg_x, deg_y, elements_x-deg_x, elements_y-deg_y, 1, &kvec_x[0], &kvec_y[0]);
  
  lrs.refine(XFIXED,0.56,0.0,0.75);
  
  std::string prefix(argv[1]);
  
  writePostscriptSuppLR(lrs, prefix);
  //writePostscriptMeshOverload(lrs,ofs);

  /*std::vector<LRBSpline2D*> bb = LinDepUtils::unpeelableBasisFunctions ( lrs );
  cout << "Is peelable: " << LinDepUtils::isPeelable(lrs) << endl;
  cout << "No. unpeelable basis functions: " << bb.size() << endl;
  return 0;
  lrs.write(cout);*/

  return 0;
}

