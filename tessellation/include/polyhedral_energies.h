#ifndef _POLYHEDRAL_ENERGIES_H
#define _POLYHEDRAL_ENERGIES_H

#include <utility>
#include <functional>
#include "common_defs.h"
#include "distance_function.h"
#include "clip_grid.h"


namespace TesselateUtils {

// ============================================================================
ValAndDer<Point2D> polygon_energy(const Point2D* const bpoints,
                                  const unsigned int num_bpoints,
                                  const Point2D* const ipoints,
                                  const unsigned int num_ipoints,
                                  const double vdist,
                                  const ClippedGrid<2>* const = nullptr); // optional
// ============================================================================


// ============================================================================
ValAndDer<Point3D> polyhedron_energy(const Point3D* const bpoints,
                                     const unsigned int num_bpoints,
                                     const Triangle* const btris,
                                     const unsigned int num_btris,
                                     const Point3D* const ipoints,
                                     const unsigned int num_ipoints,
                                     const double vdist,
                                     const ClippedGrid<3>* const = nullptr); // optional
// ============================================================================  
  

}; // end namespace Go

#endif
