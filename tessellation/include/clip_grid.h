#ifndef _CLIP_GRID_H
#define _CLIP_GRID_H

#include <array>
#include <vector>
#include "common_defs.h"

namespace TesselateUtils {

enum ClippedDomainType {
  FAR_INSIDE   = 0,    // inside polygon, 'far' from boundary
  CLOSE_INSIDE = 1,    // inside polygon, within 'vdist' from boundary
  INTERSECTED  = 2,    // polygon intersects this cell
  OUTSIDE      = 3,    // outside polygon
  UNDETERMINED = 4,    // not yet determined
};

template<int Dim>
struct ClippedGrid {
  std::array<double, 2*Dim> bbox; // bounding box [xmin, xmax, ymin, ymax, ...]
  std::array<uint, Dim> res; // resolution (number of subdivisions in each direction)
  std::array<double, Dim> cell_len; // redundant, but prevents having to recompute it all the time
  std::vector<ClippedDomainType> type; // type of each cell (is it clearly
                                       // inside or outside the polygon/shell,
                                       // or is it intersected by it).
};

// ----------------------------------------------------------------------------
ClippedGrid<2> clip_grid_polygon_2D(const Point2D* const pcorners,
                                    const uint num_pcorners,
                                    const double vdist, 
                                    const uint res_x,
                                    const uint res_y);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
ClippedGrid<3> clip_grid_shell_3D(const Point3D* const pcorners,
                                  const uint num_pcorners,
                                  const Triangle* const tris,
                                  const uint num_tris,
                                  const double vdist,
                                  const uint res_x,
                                  const uint res_y,
                                  const uint res_z);
// ----------------------------------------------------------------------------
  


  
}; // end namespace TesselateUtils

#endif
