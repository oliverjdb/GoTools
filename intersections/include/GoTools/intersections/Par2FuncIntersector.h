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

#ifndef _PAR2FUNCINTERSECTOR_H
#define _PAR2FUNCINTERSECTOR_H


#include "GoTools/intersections/IntersectorFuncConst.h"
#include "GoTools/intersections/IntersectionPoint.h"


namespace Go {


/// This class is performing intersections between a 1-dimensional
/// parametric surface and a constant.

class Par2FuncIntersector : public IntersectorFuncConst {
public:

//     Par2FuncIntersector(shared_ptr<ParamFunctionInt> func,
// 			shared_ptr<ParamFunctionInt> C,
// 			double epsge,
// 			Intersector* prev = 0);

    /// Constructor.
    /// One of the objects should refer to a 1D-surface, the other a
    /// constant (this is not checked compile-time, so we rely on the
    /// user to obey this rule).  The last two variables are relevant
    /// only if the parent has one more parameter than the Intersector
    /// to be constructed.
    /// \param func of type Param2FunctionInt.
    /// \param C of type Param0FunctionInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    Par2FuncIntersector(shared_ptr<ParamFunctionInt> func,
			shared_ptr<ParamFunctionInt> C,
			shared_ptr<GeoTol> epsge,
			Intersector *prev = 0,
			int eliminated_parameter = -1,
			double eliminated_value = 0);

    /// Destructor.
    virtual ~Par2FuncIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 2; }
	
protected:
    // Data members

    virtual shared_ptr<Intersector> 
    lowerOrderIntersector(shared_ptr<ParamFunctionInt> obj1,
			  shared_ptr<ParamFunctionInt> obj2,
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int checkCoincidence();

    virtual void microCase();
    
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    bool
    isConnected(std::vector<shared_ptr<IntersectionPoint> > bd_ints,
		int nmbbd);

    bool isConnected(std::vector<std::
		     pair<shared_ptr<IntersectionPoint>,
		     IntPtClassification> >& bd_ints, 
		     int nmb_nottouch);

    bool connectDirected(std::vector<std::
			 pair<shared_ptr<IntersectionPoint>,
			 IntPtClassification> >& bd_ints,
			 int nmbbd);

    bool canConnect(shared_ptr<IntersectionPoint> pt1,
		    shared_ptr<IntersectionPoint> pt2);

    virtual int doSubdivide();

private:

    int getSubdivisionParameter(int dir, double& par);

    // We need to decide in which direction to subdivide.
    int sortParameterDirections(int perm[]); //, int deg_edge[]);

    int checkSubdivParam(int dir, double par, double ta, double tb,
			 std::vector<shared_ptr<IntersectionPoint> >& int_pts);

    int checkIsoCurve(int pdir, bool first, double par,
 		      std::vector<shared_ptr<IntersectionPoint> > int_pts);

    bool getSubdivAtSing(int dir, double ta, double tb, double& par);

//     void splitIntResults(std::vector<std::
// 			 shared_ptr<IntersectionPoint> >& int_pts,
// 			 int pardir, double par,
// 			 double start, double end);

//     void doIterate(int pardir, double parval, double param[], double& dist,
// 		   double seed[]);

//     // Utility function for sorting input bd_int's.
//     IntPtClassification bdDir(const IntersectionPoint& int_pt,
// 			      Point sorting_dir);

    void writeDebugConnect(std::vector<std::
			   pair<shared_ptr<IntersectionPoint>,
			   IntPtClassification> >& bd_ints);

};


} // namespace Go


#endif // _PAR2FUNCINTERSECTOR_H

