#ifndef _SIMPLE_POLYHEDRON_TESSELATION_H
#define _SIMPLE_POLYHEDRON_TESSELATION_H
#include <iosfwd>
#include "TesselableVolume.h"
#include "common_defs.h"

namespace TesselateUtils
{
struct FaceLoop {
  // face orientation given by the orientation of the first edge pointed to in 'edge_ixs'
  std::vector<uint> edge_ixs; // indices to the segments making up the planar face
  bool ccw; // true if orientated counterclockwise
};

struct SimpleVolumeType {}; // we only need the type, no data necessary
  
struct SimpleSpaceTraits {
  typedef Point3D PointType;
  typedef Segment EdgeType;
  typedef FaceLoop FaceType; // indices to edges bounding the face
  typedef SimpleVolumeType VolumeType;
};

  
using SimplePolyhedron = TesselableVolume<SimpleSpaceTraits>;
  
// ----------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& os, FaceLoop fl)
// ----------------------------------------------------------------------------
{
  for (const auto e : fl.edge_ixs) os << e << ' '; os << fl.ccw << '\n';
  //for (const auto e : fl.edge_ixs) os << e << '\n';
  return os;
}

  
};

#endif
