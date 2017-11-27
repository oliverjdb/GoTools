#ifndef _GOPARAMETRICTESSELABLEVOLUME_H
#define _GOPARAMETRICTESSELABLEVOLUME_H

//#include <iostream>  
#include <tuple>
#include "common_defs.h"
#include "TesselableVolume.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/utils/Point.h"
#include "GoTools/trivariatemodel/ftVolume.h"

namespace TesselateUtils
{
  
struct GoParametricSpaceTraits {
  typedef std::shared_ptr<const Go::ParamCurve>   EdgeCurve;
  typedef std::shared_ptr<const Go::ParamSurface> FaceSurface;

  struct PointType {
    Go::Point pos; // Describes the spatial position of the point.
    Point3D param; // Describes the parameter(s) of the point relative to the
                   // unique entity that defines it (a parametric curve, surface
                   // or volume).  A given entity only defines its interior
                   // points, so e.g. boundary points of a surface are instead
                   // defined in terms of the bounding curve and its
                   // parametrization.  Corner points have no associated parameters.

    PointType() {};
    
    // construtor taking only a 3D point and setting parametrization to undefined
    PointType(const Go::Point& p)
      : pos(p), param( {std::nan(""), std::nan(""), std::nan("")} ) {};

    // constructor taking both 3D point and one parameter (useful for curves)
    PointType(const Go::Point& p, const double par)
      : pos(p), param( { par, std::nan(""), std::nan("") } ) {};

    // constructor taking 3D point and two parameters (useful for surfaces)
    PointType(const Go::Point& p, const double par1, const double par2)
      : pos(p), param( { par1, par2, std::nan("") } ) {};

    // constructor taking 3D point and three parameters (useful for volumes)
    PointType(const Go::Point& p, const Point3D& par) : pos(p), param(par) {};


  };
  
  // geometry and start/end vertex indices
  typedef std::tuple<EdgeCurve, uint, uint> EdgeType; 

  // geometry, face edge index, and forward/reverse flag
  typedef struct {
    FaceSurface surf;
    vector<pair<uint, bool>> ix; // boolean indicate orientation of edge.  If
                                 // 'false', edge is oriented in the reverse
                                 // direction
  } FaceType;

  //typedef std::pair<FaceSurface, vector<std::pair<uint, bool>>> FaceType; 
  typedef std::shared_ptr<const Go::ParamVolume>  VolumeType;
};

  
using GoParametricTesselableVolume = TesselableVolume<GoParametricSpaceTraits>;  

// ----------------------------------------------------------------------------
// Additional constructor, making a TesselableVolume out of a Go::ftVolume
template<> template<>
GoParametricTesselableVolume::TesselableVolume(Go::ftVolume& fvol);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& os,
                                const GoParametricSpaceTraits::PointType& p)
// ----------------------------------------------------------------------------
{
  // writing only 3D part (not parameters) for now
  p.pos.write(os);

  // if you want to write the parameterization, uncomment the following line
  //os << p.param[0] << " " << p.param[1] << " " << p.param[2] << '\n';
  return os;
}


  
};

#endif
