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

#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <iostream>
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 4) {
	//MESSAGE("Usage: input_pts tolerance output_g2 spline_deg");
	MESSAGE("Usage: input_pts tolerance output_g2");
	return 0;
    }

    // Read input arguments
    std::ifstream inpts(argv[1]);
    double epsge = atof(argv[2]);
    std::ofstream outg2(argv[3]);
    //int deg = atoi(argv[4]);

    ALWAYS_ERROR_IF(inpts.bad(), "Input file not found or file corrupt");
    
    vector<double> pts,pars;
    double tmpx,tmpy;

    int count = 0;
    while(!inpts.eof()) {
        inpts >> tmpx;
	inpts >> tmpy;
	if (inpts.eof()) break;
        //cout << tmpx << " " << tmpy << endl;
	pts.push_back(tmpx);
	pts.push_back(tmpy);
	count += 1;
        //pars.push_back((double)count); // for uniform parametrization
    }

    // for chordal parametrization
    pars.push_back(0.0);
    for (int ix=1; ix!=count; ++ix) {
        pars.push_back(pars[ix-1] + sqrt(sqrt(pow((pts[2*(ix-1)]-pts[2*ix]),2) + pow((pts[2*(ix-1)+1]-pts[2*ix+1]),2))));
    }

    //for (int ix=0; ix!=count; ++ix) cout << pts[2*ix] << " " << pts[2*ix+1] << " " << pars[ix] << endl;

    const int dim = 2;
    ApproxCurve appcrv(pts,pars,dim,epsge);

    //cout << "Tolerance: " << epsge << endl;
    //cout << "Dimension: " << dim << endl;

    double maxdist, avdist, max_iter = 10;

    vector<Point> start = {Point(pts[0],pts[1])};
    vector<Point> end = {Point(pts[pts.size()-2],pts[pts.size()-1])};

    appcrv.setEndPoints(start,end);
    shared_ptr<SplineCurve> curve = appcrv.getApproxCurve(maxdist, avdist, max_iter);
   
    //cout << "Max dist: " << maxdist << endl;
    //cout << "Av. dist: " << avdist << endl;

    //curve->writeStandardHeader(cout);
    //curve->write(cout);

    cout << count << ", " << curve->numCoefs() << ", " << maxdist << ", " << avdist << endl;

    outg2 << "{\"control_pts\": [";
    auto it=curve->coefs_begin();
    for (; it!=curve->coefs_end()-2; it+=2) {
        outg2 << "[" << (int) it[0] << ", " << (int) it[1] << "], "; 
    }
    outg2 << "[" << (int) it[0] << ", " << (int) it[1] << "]],\n";
    outg2 << "\"knots\": [";
    auto itk=curve->knotsBegin();
    for (; itk!=curve->knotsEnd()-1; ++itk) {
        outg2 << *itk << ", ";
    }
    outg2 << *itk << "]\n}";

}

