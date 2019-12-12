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

#include "GoTools/creators/SmoothSurf.h"
#include <fstream>
#include <limits>     

using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 5) {
	MESSAGE("Usage: ./makeSplineApproximations inputfile outputfile numcoefs_x numcoefs_y");
	return 0;
    }

    // Read input arguments
    std::ifstream ifs(argv[1]);
    ALWAYS_ERROR_IF(ifs.bad(), "Input file not found or file corrupt");

    std::ofstream ofs(argv[2]);

    vector<double> param;
    vector<double> data;

    double x,y,h;
    // Input surface should be a list of points, x,y,h
    double xmin = numeric_limits<double>::max();
    double xmax = numeric_limits<double>::min();
    double ymin = numeric_limits<double>::max();
    double ymax = numeric_limits<double>::min();

    ifs >> x >> y >> h;
    while (true) {
	if (x < xmin) xmin = x;
	if (x > xmax) xmax = x;
	if (y < ymin) ymin = y;
	if (y > ymin) ymax = y;
	param.push_back(x);
	param.push_back(y);
	//data.push_back(x); // uncomment if wanting 3D spline
	//data.push_back(y); 
	data.push_back(h);
        ifs >> x >> y >> h;
	if (ifs.eof()) break;
    }


    // Hard coded 14x14
    int num_coefs_x = atoi(argv[3]);  
    int num_coefs_y = atoi(argv[4]); 
    
    int deg_x = 2;
    int deg_y = 2;
    int dim = 1;
    
    double xstep = (xmax - xmin)/(num_coefs_x-deg_x);
    double ystep = (ymax - ymin)/(num_coefs_y-deg_y);

    vector<double> knots_x(num_coefs_x+deg_x+1, 0.0);
    vector<double> knots_y(num_coefs_y+deg_y+1, 0.0);
    for (int ix=0; ix!=deg_x+1; ++ix) {
        knots_x[ix] = xmin;
	knots_x[num_coefs_x+deg_x-ix] = xmax;
    }
    for (int ix=0; ix!=deg_y+1; ++ix) {
        knots_y[ix] = ymin;
        knots_y[num_coefs_y+deg_y-ix] = ymax;
    }
    for (int ix=1; ix!=num_coefs_x-deg_x; ++ix) {
        knots_x[deg_x+ix] = xmin + xstep*ix;
    }
    for (int ix=1; ix!=num_coefs_y-deg_y; ++ix) {
        knots_y[deg_y+ix] = ymin + ystep*ix;
    }

    //for (int ix=0; ix!=num_coefs_x+deg_x+1; ++ix) {
    //    cout << knots_x[ix] << " ";
    //} cout << endl;
    //for (int ix=0; ix!=num_coefs_y+deg_y+1; ++ix) {
    //    cout << knots_y[ix] << " ";
    //} cout << endl;

    vector<int> zeros(num_coefs_x*num_coefs_y, 0);
    vector<double> weights(data.size());
    vector<double> ones(data.size(), 1.0);
    for (int ix=0; ix!=data.size(); ++ix) {
        weights[ix] = pow(1.4,-1.0*fabs(data[ix])); // most weight at distance 0, tails off exponentially
	//cout << data[ix] << " " << weights[ix] << endl; 
    }
    
    shared_ptr<SplineSurface> surf(new SplineSurface(num_coefs_x, num_coefs_y, deg_x+1, deg_y+1, knots_x.begin(), knots_y.begin(), zeros.begin(), dim));
    
    SmoothSurf approx;
    approx.attach(surf, &zeros[0], &zeros[0]);
    approx.setLeastSquares(data, param, weights, 0.999);
    approx.equationSolve(surf);

    int count = 0;
    for (auto it=surf->coefs_begin(); it!=surf->coefs_end(); it++) {
        count++;
	ofs << (*it) << " ";
        if (count % num_coefs_x == 0) ofs << endl;	
    } 

}

