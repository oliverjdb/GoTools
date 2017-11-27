#ifndef _TESSELABLE_VOLUME_H
#define _TESSELABLE_VOLUME_H

#include <iosfwd>
#include <vector>
#include <functional>// @@ only for marking unimplemented functions
#include <stdexcept> // @@ only for marking unimplemented functions
#include <iostream>  // @@ only for debug purposes
#include <iterator>
#include "common_defs.h"
#include "tesselate_utils.h"

namespace TesselateUtils {

// ============================================================================
// The template type should provide the four typedefs:
// PointType, EdgeType, FaceType and VolumeType  
template<typename BoundedSpaceTraits>
class TesselableVolume
// ============================================================================
{
public:
  typedef typename BoundedSpaceTraits::PointType  PointType;
  typedef typename BoundedSpaceTraits::EdgeType   EdgeType;
  typedef typename BoundedSpaceTraits::FaceType   FaceType;
  typedef typename BoundedSpaceTraits::VolumeType VolumeType;
  
  TesselableVolume(const std::vector<PointType>& c,
                     const std::vector<EdgeType>& e,
                     const std::vector<FaceType>& f,
                     const VolumeType& v) :
    corners_(c), edges_(e), faces_(f), volume_(v)
  {}

  // Possible alternative constructor
  template<typename ArgumentClass> TesselableVolume(ArgumentClass& arg);
      
  std::ostream& write (std::ostream& os) const;
  void writeTesselatedOutline(std::ostream& os) const;
  void writeTesselatedShell(std::ostream& os) const;
  void writeTesselatedVolume(std::ostream& os) const;
  
  std::vector<PointType> cornerPoints() const;
  std::vector<PointType> edgePoints() const; // all edge points
  std::vector<PointType> edgePoints(uint edge_ix) const;
  std::vector<PointType> facePoints() const;
  std::vector<PointType> facePoints(uint face_ix) const;
  std::vector<PointType> volumePoints() const;

  std::vector<Tet> getTets() const { return volume_tets_;}
  
  void tesselate(const double vdist); 
  bool is_tesselated() const {return !edge_ipoints_.empty();}
  
  uint numCorners() const { return (uint)corners_.size();}
  uint numEdges()   const { return (uint)edges_.size();}
  uint numFaces()   const { return (uint)faces_.size();}

  uint numEdgePoints() const;
  uint numFacePoints() const;
  uint numVolumePoints() const;

  uint numTets() const { return (uint)volume_tets_.size();}

  std::array<PointType, 2> edgeCorners(uint edge_ix) const;
  std::vector<PointType> faceBoundaryPoints(uint face_ix) const;
  
  // the following functions need template specialization
  std::array<uint, 2> edgeCornerIndices(uint edge_ix) const;
  std::vector<uint> faceBoundaryPointIndices(uint face_ix) const;
    
private:
  void update_triangle_indices();
  void update_tet_indices();
  void compute_global_edge_point_indices();
  void compute_global_face_point_indices();

  // The following tesselating functions must be specifically implemented
  // for each particular choice of BoundedSpaceTraits.
  static void compute_tesselation(const std::array<PointType, 2>& boundary,
                                  const EdgeType& edge,
                                  const double vdist,
                                  std::vector<PointType>& ipoints);
  static void compute_tesselation(const std::vector<PointType>& boundary,
                                  const FaceType& face,
                                  const double vdist,
                                  std::vector<PointType>& ipoints,
                                  std::vector<Triangle>& triangles);
  static void compute_tesselation(const std::vector<PointType>& bpoints,
                                  const std::vector<Triangle>& btris,
                                  const VolumeType& volume,
                                  const double vdist,
                                  std::vector<PointType>& ipoints,
                                  std::vector<Tet>& tets);
  
  // the following properties always have to be set.  They define the boundary
  // of the volume
  std::vector<PointType> corners_;
  std::vector<EdgeType>  edges_;
  std::vector<FaceType>  faces_;
  VolumeType volume_;

  // the following propeties pertains to the tesselation, and will be created
  // upon request
  std::vector<std::vector<PointType>> edge_ipoints_;
  std::vector<std::vector<PointType>> face_ipoints_;
  std::vector<PointType>              volume_ipoints_;
  std::vector<std::vector<Triangle>>  face_triangles_;
  std::vector<Tet>                    volume_tets_;

  std::vector<uint> edge_ipoints_start_ixs_;  // global start index of the
                                              // corresponding edge's internal
                                              // points
  std::vector<uint> face_ipoints_start_ixs_;  // global start index of the
                                              // corresponding face's internal
                                              // points
}; // end struct TesselableVolume


// ----------------------------------------------------------------------------        
template<typename BoundedSpaceTraits> inline void
TesselableVolume<BoundedSpaceTraits>::compute_global_edge_point_indices() 
// ----------------------------------------------------------------------------
{
  edge_ipoints_start_ixs_.resize(numEdges(), numCorners());
  for (uint i = 0; i != numEdges()-1; ++i)
    edge_ipoints_start_ixs_[i+1] = edge_ipoints_start_ixs_[i] +
                                   (uint)edge_ipoints_[i].size();
}

// ----------------------------------------------------------------------------        
template<typename BoundedSpaceTraits> inline void
TesselableVolume<BoundedSpaceTraits>::compute_global_face_point_indices() 
// ----------------------------------------------------------------------------
{
  face_ipoints_start_ixs_.resize(numFaces(), (uint)edgePoints().size());
  for (uint i = 0; i != numFaces() - 1; ++i)
    face_ipoints_start_ixs_[i+1] = face_ipoints_start_ixs_[i] +
                                   (uint)face_ipoints_[i].size();
}

  
// ----------------------------------------------------------------------------        
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::faceBoundaryPoints(uint face_ix) const
// ----------------------------------------------------------------------------
{
  const std::vector<PointType> pts = edgePoints();
  const std::vector<uint> ixs = faceBoundaryPointIndices(face_ix);
  std::vector<PointType> result(ixs.size());
  transform(ixs.begin(), ixs.end(), result.begin(),
            [&pts] (uint ix) { return pts[ix];});
  return result;
}
  
// ----------------------------------------------------------------------------      
template<typename BoundedSpaceTraits> inline
uint TesselableVolume<BoundedSpaceTraits>::numEdgePoints() const
// ----------------------------------------------------------------------------
{
  return numCorners() +
    accumulate(edge_ipoints_.begin(), edge_ipoints_.end(), 0,
               [](uint acc, const std::vector<PointType>& v)
               { return acc + (uint)v.size();});
}

// ----------------------------------------------------------------------------      
template<typename BoundedSpaceTraits> inline
uint TesselableVolume<BoundedSpaceTraits>::numFacePoints() const
// ----------------------------------------------------------------------------
{
  return numEdgePoints() +
    accumulate(face_ipoints_.begin(), face_ipoints_.end(), 0, 
               [](uint acc, const std::vector<PointType>& v)
               { return acc + (uint)v.size();});
}

// ----------------------------------------------------------------------------      
template<typename BoundedSpaceTraits> inline
uint TesselableVolume<BoundedSpaceTraits>::numVolumePoints() const
// ----------------------------------------------------------------------------
{
  return numFacePoints() + (uint)volume_ipoints_.size();
}
  
// ----------------------------------------------------------------------------    
template<typename BoundedSpaceTraits> inline
void TesselableVolume<BoundedSpaceTraits>::update_triangle_indices()
// ----------------------------------------------------------------------------    
{
  for (uint f_ix = 0; f_ix != numFaces(); ++f_ix) {
    
    const auto glob_bpoint_indices = faceBoundaryPointIndices(f_ix);
    const uint num_bpoints = (uint)glob_bpoint_indices.size();
    const uint ipoints_ix_start = face_ipoints_start_ixs_[f_ix];

    for (auto& t : face_triangles_[f_ix])
      // looping over triangle corners and updating indices
      for (uint i = 0; i != 3; ++i) 
        t[i] = ( t[i] < num_bpoints )  ?
               glob_bpoint_indices[t[i]] : 
               (t[i] - num_bpoints) + ipoints_ix_start;
  }
}

// // ----------------------------------------------------------------------------    
// template<typename BoundedSpaceTraits> inline
// void TesselableVolume<BoundedSpaceTraits>::update_tet_indices()
// // ----------------------------------------------------------------------------    
// {
//   std::cout << "Warning: update_tet_indices() not yet implemented.  Skipping." << std::endl;
// }
  
// ----------------------------------------------------------------------------  
template<typename BoundedSpaceTraits>
inline void TesselableVolume<BoundedSpaceTraits>::tesselate(const double vdist)
// ----------------------------------------------------------------------------
{
  if (is_tesselated())
    return; // nothing more to do

  edge_ipoints_.resize(numEdges());
  face_ipoints_.resize(numFaces());
  face_triangles_.resize(numFaces());
  //volume_tets_.resize(numFaces());

  // tesselate all edges
  for (uint i = 0; i != numEdges(); ++i)
    compute_tesselation(edgeCorners(i), edges_[i], vdist, edge_ipoints_[i]);
  compute_global_edge_point_indices();
  
  // tesselate all faces
  for (uint i = 0; i != numFaces(); ++i)
    compute_tesselation(faceBoundaryPoints(i), faces_[i], vdist,
                        face_ipoints_[i], face_triangles_[i]);
  compute_global_face_point_indices();
  

  // make globally consistent indexing of points for the triangles
  update_triangle_indices();

  // tesselate the volume
  compute_tesselation(facePoints(), flatten(face_triangles_),
                      volume_, vdist, volume_ipoints_, volume_tets_);

  // remove already-existing boundary points from volume_ipoints_
  volume_ipoints_ = std::vector<PointType>(volume_ipoints_.begin() + numFacePoints(),
                                           volume_ipoints_.end());
  
  // // make globally consistent indexing of points for the triangles
  // update_tet_indices();
}
  
// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits>  
std::ostream& TesselableVolume<BoundedSpaceTraits>::write(std::ostream& os) const
// ----------------------------------------------------------------------------
{
  // streaming basic information
  for (const auto c : corners_) os << c;  os << '\n';
  for (const auto e : edges_)   os << e;  os << '\n';
  for (const auto f : faces_)   os << f;  os << '\n';

  // streaming tesselation information, if any
  for (const auto eip : edge_ipoints_)
    for (const auto ip : eip) os << ip; os << '\n';
  for (const auto fip : face_ipoints_)
    for (const auto ip : fip) os << ip; os << '\n';
  for (const auto ip : volume_ipoints_) os << ip; os << '\n';
  for (const auto ft : face_triangles_)
    for (const auto t : ft) os << t; os << '\n';
  for (const auto vt : volume_tets_)
    for (const auto t : vt) os << t; os << '\n';
  return os;
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
std::ostream&
operator<<(std::ostream& os, const TesselableVolume<BoundedSpaceTraits>& tv)
// ----------------------------------------------------------------------------
{
  return tv.write(os);
}

// ----------------------------------------------------------------------------  
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::cornerPoints() const
// ----------------------------------------------------------------------------  
{
  return corners_;
}

// ----------------------------------------------------------------------------  
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::edgePoints() const
// ----------------------------------------------------------------------------  
{
  std::vector<PointType> tmp(corners_);
  for (const auto& e : edge_ipoints_)
    tmp.insert(tmp.end(), e.begin(), e.end());
  return tmp;
}

// ----------------------------------------------------------------------------  
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::facePoints() const
// ----------------------------------------------------------------------------  
{
  std::vector<PointType> tmp(edgePoints());
  for (const auto& f : face_ipoints_)
    tmp.insert(tmp.end(), f.begin(), f.end());
  return tmp;
}

// ----------------------------------------------------------------------------  
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::volumePoints() const
// ----------------------------------------------------------------------------
{
  std::vector<PointType> tmp(facePoints());
  tmp.insert(tmp.end(), volume_ipoints_.begin(), volume_ipoints_.end());
  return tmp;
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::edgePoints(uint edge_ix) const
// ----------------------------------------------------------------------------
{
  const auto endpoints = edgeCorners(edge_ix);
  std::vector<PointType> result {endpoints.first, endpoints.second};
  result.insert(result.end(),
                edge_ipoints_[edge_ix].begin(),
                edge_ipoints_[edge_ix].end());
  return result;
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
std::vector<typename TesselableVolume<BoundedSpaceTraits>::PointType>
TesselableVolume<BoundedSpaceTraits>::facePoints(uint face_ix) const
// ----------------------------------------------------------------------------  
{
  std::vector<PointType> result(faceBoundaryPoints(face_ix));
  result.insert(result.end(),
                face_ipoints_[face_ix].begin(),
                face_ipoints_[face_ix].end());
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
std::array<typename TesselableVolume<BoundedSpaceTraits>::PointType, 2>
TesselableVolume<BoundedSpaceTraits>::edgeCorners(uint edge_ix) const
// ----------------------------------------------------------------------------
{
  const auto ixs = edgeCornerIndices(edge_ix);
  return {corners_[ixs[0]], corners_[ixs[1]]};
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
void TesselableVolume<BoundedSpaceTraits>::writeTesselatedOutline(std::ostream& os) const
// ----------------------------------------------------------------------------  
{
  // write all tesselated edges on a easily plottable format

  // first, write all corner points
  for (const auto& c: corners_) os << c << '\n';

  // then write all the interior points
  for (const auto& pvec : edge_ipoints_)
    for (const auto& p : pvec) os << p << '\n';
  os << '\n';
  
  // then write the connectivities
  for (uint edge_ix = 0; edge_ix != numEdges(); ++edge_ix) {

    // write first start point index
    os << edgeCornerIndices(edge_ix)[0] << ' ';

    // write end point index, and start point index of next segment (which is
    // the same point)
    if (is_tesselated()) {
      const uint start_ix = edge_ipoints_start_ixs_[edge_ix];
      for (uint ip_ix = 0; ip_ix != (uint)edge_ipoints_[edge_ix].size(); ++ip_ix)
        os << start_ix + ip_ix << '\n' << start_ix + ip_ix << ' ';
    }

    // write last end point index
    os << edgeCornerIndices(edge_ix)[1] << '\n';
    
  }
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
void TesselableVolume<BoundedSpaceTraits>::writeTesselatedShell(std::ostream& os) const
// ----------------------------------------------------------------------------
{
  // writing points
  os << "Points: " << '\n';
  const auto shell_points = facePoints();
  std::copy(shell_points.begin(), shell_points.end(), std::ostream_iterator<PointType>(os, " "));

  
  // writing triangles
  os << "\n" << "Triangles: " << '\n';
  for (const auto& tris : face_triangles_)
    for (const auto& t : tris)
      os << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
}

// ----------------------------------------------------------------------------
template<typename BoundedSpaceTraits> inline
void TesselableVolume<BoundedSpaceTraits>::writeTesselatedVolume(std::ostream& os) const
// ----------------------------------------------------------------------------
{
  // writing points
  os << "Points: " << '\n';
  const auto vpoints  = volumePoints();
  std::copy(vpoints.begin(), vpoints.end(), std::ostream_iterator<PointType>(os, " "));

  os << std::endl;;
  // writing tets
  os << "Tets: " << '\n';
  for (const auto& t : volume_tets_)
    os << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << '\n';
}



}; // end namespace TesselateUtils


#endif

