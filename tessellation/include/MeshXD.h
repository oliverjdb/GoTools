#ifndef _MESH_XD_H
#define _MESH_XD_H

#include "common_defs.h"

namespace TesselateUtils
{
  struct Mesh2D {
    std::vector<Point2D> points;
    std::vector<Triangle> tris;
  };

  struct Mesh3D {
    std::vector<Point3D> points;
    std::vector<Tet> tets;
  };
}

#endif
