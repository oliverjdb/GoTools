#ifndef _COMMON_DEFS_H
#define _COMMON_DEFS_H

#include <array>
#include <vector>
#include <cmath>
#include <ostream>
// #include "GoTools/utils/Point.h"

namespace TesselateUtils {

using uint	= unsigned int;
using Segment   = std::array<uint, 2>; // Representing a segment by indices to endpoints
using Triangle	= std::array<uint, 3>;
using Tet	= std::array<uint, 4>;

template<int Dim>
struct PointXD {

  double& operator[](uint i) { return coords[i];}
  const double& operator[](uint i) const {return coords[i];}
  template<typename Iterator>
  void copyFrom(Iterator begin) {std::copy(begin, begin+Dim, &coords[0]);}
  
  std::array<double, Dim> coords;
};

using Point1D = PointXD<1>;
using Point2D = PointXD<2>; //std::array<double, 2>;
using Point3D = PointXD<3>; // std::array<double, 3>;

template<typename PointXD>  
struct ValAndDer {
  double val;
  std::vector<PointXD> der;

  void reset(uint num_der);
};

template<> inline void ValAndDer<Point2D>::reset(uint num_der) {
  val = 0;
  der = std::vector<Point2D>(num_der, Point2D {0, 0});
}

  template<> inline void ValAndDer<Point3D>::reset(uint num_der) {
  val = 0;
  der = std::vector<Point3D>(num_der, Point3D {0, 0, 0});
}

// // ========================= Adapter functions for Go =========================

// // ----------------------------------------------------------------------------
// inline double dist2(const Go::Point& p1, const Go::Point& p2)
// // ----------------------------------------------------------------------------
// {
//   return p1.dist2(p2);
// }
  
// =============================== 1D operators ===============================

// ----------------------------------------------------------------------------
inline void operator *= (Point1D& p, const double d)
// ----------------------------------------------------------------------------
{
  p[0] *= d;
}
  
// ----------------------------------------------------------------------------
inline Point1D operator*(const Point1D& p, double t)
// ----------------------------------------------------------------------------
{
  Point1D tmp(p);
  tmp *= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point1D operator*(double t, const Point1D& p)
// ----------------------------------------------------------------------------
{
  Point1D tmp(p);
  tmp *= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline void operator += (Point1D& p1, const Point1D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] += p2[0]; 
}

// ----------------------------------------------------------------------------
inline void operator -= (Point1D& p1, const Point1D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] -= p2[0]; 
}

  

// =============================== 2D operators ===============================
  
// ----------------------------------------------------------------------------
inline void operator += (Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] += p2[0]; 
  p1[1] += p2[1];
}

// ----------------------------------------------------------------------------
inline void operator -= (Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] -= p2[0]; 
  p1[1] -= p2[1];
}

// ----------------------------------------------------------------------------
inline void operator *= (Point2D& p1, const double t)
// ----------------------------------------------------------------------------
{ 
  p1[0] *= t; 
  p1[1] *= t;
}

// ----------------------------------------------------------------------------
inline void operator /= (Point2D& p1, const double t)
// ----------------------------------------------------------------------------
{ 
  p1[0] /= t; 
  p1[1] /= t;
}

// ----------------------------------------------------------------------------
inline Point2D operator + (const Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p1);
  tmp += p2;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point2D operator - (const Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p1);
  tmp -= p2;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point2D operator*(const Point2D& p, double t)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p);
  tmp *= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point2D operator/(const Point2D& p, double t)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p);
  tmp /= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& os, const Point2D& p)
// ----------------------------------------------------------------------------
{
  os << p[0] << ' ' << p[1] << '\n';
  return os;
}

// ----------------------------------------------------------------------------
inline double norm2(const Point2D& p)
// ----------------------------------------------------------------------------
{
  return p[0] * p[0] + p[1] * p[1];
}

// ----------------------------------------------------------------------------
  inline double norm(const Point2D& p)
// ----------------------------------------------------------------------------
{
  return std::sqrt(norm2(p));
}
  
// ----------------------------------------------------------------------------
inline double dist2(const Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{
  const double dx = p1[0] - p2[0];
  const double dy = p1[1] - p2[1];
  return dx*dx + dy*dy;
}

// ----------------------------------------------------------------------------
inline bool acute_angle(const Point2D& a, const Point2D& b, const Point2D& c)
// ----------------------------------------------------------------------------
{
  // angle is acute if scalar product of vectors ba and bc is positive
  return (c[0] - b[0]) * (a[0] - b[0]) + (c[1] - b[1]) * (a[1] - b[1]) > 0;
}

// ----------------------------------------------------------------------------
// scalar product
inline double operator*(const Point2D& p, const Point2D& q)
// ----------------------------------------------------------------------------
{
  return p[0] * q[0] + p[1] * q[1];
}

  
// =============================== 3D operators ===============================

// ----------------------------------------------------------------------------
inline void operator += (Point3D& p1, const Point3D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] += p2[0]; 
  p1[1] += p2[1];
  p1[2] += p2[2];
}

// ----------------------------------------------------------------------------
inline void operator -= (Point3D& p1, const Point3D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] -= p2[0]; 
  p1[1] -= p2[1];
  p1[2] -= p2[2];
}

// ----------------------------------------------------------------------------
inline void operator *= (Point3D& p1, const double t)
// ----------------------------------------------------------------------------
{ 
  p1[0] *= t; 
  p1[1] *= t;
  p1[2] *= t;
}

// ----------------------------------------------------------------------------
inline void operator /= (Point3D& p1, const double t)
// ----------------------------------------------------------------------------
{ 
  p1[0] /= t; 
  p1[1] /= t;
  p1[2] /= t;
}

// ----------------------------------------------------------------------------
inline Point3D operator + (const Point3D& p1, const Point3D& p2)
// ----------------------------------------------------------------------------
{
  Point3D tmp(p1);
  tmp += p2;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point3D operator - (const Point3D& p1, const Point3D& p2)
// ----------------------------------------------------------------------------
{
  Point3D tmp(p1);
  tmp -= p2;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point3D operator*(const Point3D& p, double t)
// ----------------------------------------------------------------------------
{
  Point3D tmp(p);
  tmp *= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point3D operator*(double t, const Point3D& p)
// ----------------------------------------------------------------------------
{
  Point3D tmp(p);
  tmp *= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point3D operator/(const Point3D& p, double t)
// ----------------------------------------------------------------------------
{
  Point3D tmp(p);
  tmp /= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& os, const Point3D& p)
// ----------------------------------------------------------------------------
{
  os << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
  return os;
}
 
// ----------------------------------------------------------------------------
 inline double norm2(const Point3D& p)
// ----------------------------------------------------------------------------
{
  return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

// ----------------------------------------------------------------------------
  inline double norm(const Point3D& p)
// ----------------------------------------------------------------------------
{
  return std::sqrt(norm2(p));
}
  
// ----------------------------------------------------------------------------
inline double dist2(const Point3D& p1, const Point3D& p2)
// ----------------------------------------------------------------------------
{
  const double dx = p1[0] - p2[0];
  const double dy = p1[1] - p2[1];
  const double dz = p1[2] - p2[2];
  return dx*dx + dy*dy + dz*dz;
}


// ----------------------------------------------------------------------------
inline bool acute_angle(const Point3D& a, const Point3D& b, const Point3D& c)
// ----------------------------------------------------------------------------
{
  // angle is acute if scalar product of vectors ba and bc is positive
  return (c[0] - b[0]) * (a[0] - b[0]) +
         (c[1] - b[1]) * (a[1] - b[1]) +
         (c[2] - b[2]) * (a[2] - b[2])> 0;
}

// ----------------------------------------------------------------------------
// cross product
inline Point3D operator^(const Point3D& p, const Point3D& q)
// ----------------------------------------------------------------------------  
{
  return { p[1] * q[2] - p[2] * q[1],
           p[2] * q[0] - p[0] * q[2],
           p[0] * q[1] - p[1] * q[0]};
}


// ----------------------------------------------------------------------------
// scalar product
inline double operator*(const Point3D& p, const Point3D& q)
// ----------------------------------------------------------------------------
{
  return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
}

// ================= Operators on segments, triangles and tets =================

inline std::ostream& operator<<(std::ostream& os, const Segment& s)
{
  os << s[0] << ' ' << s[1] << '\n';
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
  os << tri[0] << ' ' << tri[1] << ' ' << tri[2] << '\n';
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const Tet& tet)
{
  os << tet[0] << ' ' << tet[1] << ' ' << tet[2] << ' ' << tet[3] << '\n';
  return os;
}


  
}; // end namespace TesselateUtils

#endif
