#ifndef _BASIC_INTERSECTIONS_H
#define _BASIC_INTERSECTIONS_H

#include "common_defs.h"

namespace TesselateUtils {

enum IsectCase {
  DISJOINT,         // no intersection
  AT_POINT,         // share a common point only
  ALONG_BOUNDARY,   // boundary of one intersects interior or boundary of other
  OVERLAPPING,      // interiors of A and B intersect
  CONTAINING,       // one fully contained in other
  COLINEAR_OVERLAP  // two 2D segments overlap colinearly 
};
  
// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT or CONTAINING
IsectCase isect_point_seg_1D(const double p, const double* const s, const double tol);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT or CONTAINING
IsectCase isect_point_seg_2D(const Point2D& p, const Point2D* const s, const double tol);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT or CONTAINING
IsectCase isect_point_seg_3D(const Point3D& p, const Point3D* const s, const double tol);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, CONTAINING or OVERLAPPING
IsectCase isect_seg_seg_1D(const double* const s1,
                            const double* const s2,
                            const double tol,
                            double* const result);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, COLINEAR_OVERLAP, CONTAINING or OVERLAPPING  
IsectCase isect_seg_seg_2D(const Point2D* const s1,
                           const Point2D* const s2,
                           const double tol,
                           Point2D* const result);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, COLINEAR_OVERLAP, CONTAINING or OVERLAPPING  
IsectCase isect_seg_seg_3D(const Point3D* const s1,
                           const Point3D* const s2,
                           const double tol,
                           Point3D* const result);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, ALONG_BOUDNARY, CONTAINING or OVERLAPPING    
IsectCase isect_seg_triangle_2D(const Point2D* const s, // pointer to two points
                                const Point2D* const tri, // pointer to three points
                                const double tol,
                                Point2D* const result);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
IsectCase isect_line_triangle_2D(const Point2D* const s, // segment defining the line
                                 const Point2D* const tri,
                                 const double tol,
                                 Point2D* const result);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, ALONG_BOUNDARY, CONTAINING or OVERLAPPING
IsectCase isect_triangle_triangle_2D(const Point2D* const tri1,
                                     const Point2D* const tri2,
                                     const double tol);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, ALONG_BOUNDARY, CONTAINING or OVERLAPPING
IsectCase isect_triangle_triangle_3D(const Point3D* const tri1,
                                     const Point3D* const tri2,
                                     const double tol);
// ----------------------------------------------------------------------------    
  
}; // end namespace TesselateUtils


#endif
