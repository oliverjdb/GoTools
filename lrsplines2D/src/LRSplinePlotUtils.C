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

#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/lrsplines2D/Mesh2DIterator.h"

namespace {

  void highlight_non_unitary_gamma(Go::LRSplineSurface& lr_spline_sf, std::ostream &out)
  {
    for ( Go::LRSplineSurface::BSplineMap::const_iterator it = lr_spline_sf.basisFunctionsBegin(); 
	  it != lr_spline_sf.basisFunctionsEnd(); ++it ) {
      if ( it->second->gamma() < 0.9999999 ) {
	out << "newpath\n";
	out << it->second->umin() << " " << it->second->vmin() << " moveto\n";
	out << it->second->umax() << " " << it->second->vmax() << " lineto\n";
	out << "stroke\n";
      }  
    }
    out << "%%EOF\n";
  }

}

namespace Go
{
  void writeMesh(const Go::Mesh2D& mesh, std::ostream &out, const double scale)
  {
    const bool close = false;
    const bool colorDiag = false;
    const int num_diff_knots_u = mesh.numDistinctKnots(XFIXED);
    const int num_diff_knots_v = mesh.numDistinctKnots(YFIXED);
    // std::vector<double> knot_u(num_diff_knots_u), knot_v(num_diff_knots_v);
    const double* const knot_u = mesh.knotsBegin(XFIXED);
    const double* const knot_v = mesh.knotsBegin(YFIXED);
    double min_span_u = knot_u[1] - knot_u[0];
    double min_span_v = knot_v[1] - knot_v[0];
    for(size_t i=1; i<num_diff_knots_u-1; i++)
      min_span_u = (min_span_u<knot_u[i+1]-knot_u[i]) ? min_span_u : knot_u[i+1]-knot_u[i];
    for(size_t i=1; i<num_diff_knots_v-1; i++)
      min_span_v = (min_span_v<knot_v[i+1]-knot_v[i]) ? min_span_v : knot_v[i+1]-knot_v[i];

    // get date
    time_t t = time(0);
    tm* lt = localtime(&t);
    char date[11];
    sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

    // get bounding box
    double umin = mesh.minParam(XFIXED);
    double umax = mesh.maxParam(XFIXED);
    double vmin = mesh.minParam(YFIXED);
    double vmax = mesh.maxParam(YFIXED);
    double dx = umax - umin;
    double dy = vmax - vmin;
    //double scale = 1.0;//(dx>dy) ? 1000.0/dx : 1000.0/dy;
    // set the duplicate-knot-line (dkl) display width
    double dkl_range = (min_span_u>min_span_v) ? min_span_v*scale/6.0 : min_span_u*scale/6.0; 
    double xmin = (umin - dx/100.0)*scale;
    double ymin = (vmin - dy/100.0)*scale;
    double xmax = (umax   + dx/100.0)*scale;// + dkl_range;
    double ymax = (vmax   + dy/100.0)*scale;// + dkl_range;

    // print eps header
    out << "%!PS-Adobe-3.0 EPSF-3.0\n";
    out << "%%Creator: LRSplinePlotUtils.C object\n";
    out << "%%Title: LRSplineSurface parameter domain\n";
    out << "%%CreationDate: " << date << std::endl;
    out << "%%Origin: 0 0\n";
    out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;

    out << "0 setgray\n";
    //	double linewidth = std::min(dx/100.0, dy/100.0);
    double linewidth = 0.1*scale*std::min(dx/(10.0*num_diff_knots_u), dy/(10.0*num_diff_knots_v));
    out << linewidth << " setlinewidth\n";
    Mesh2DIterator mesh_beg = mesh.begin();
    Mesh2DIterator mesh_end = mesh.end();
    Mesh2DIterator mesh_it = mesh_beg;
    while (mesh_it != mesh_end)
      {
	const std::array<int, 4> elem_corners = *mesh_it; // ll_u, ll_v, ur_u, ur_v.
	// Currently we do not care about multiplicity.
	double dm= 0.0;// (mesh_iter[i]->multiplicity_==1) ? 0 : dkl_range/(mesh_iter[i]->multiplicity_-1);
	double m = 1.0; // The multiplicity,
	//	    int mult = mesh_it->multiplicity_;
	// First we create the lines in the u-dir.
	// We also do not care about double lines (for neighbour elements) ...
	// out << mesh_iter[i]->start_*scale << " " << mesh_iter[i]->const_par_*scale + dm*m << " moveto\n";
	// if(mesh_iter[i]->stop_ == umax)
	//     out << mesh_iter[i]->stop_*scale+dkl_range << " " << mesh_iter[i]->const_par_*scale + dm*m << " lineto\n";
	// else
	//     out << mesh_iter[i]->stop_*scale << " " << mesh_iter[i]->const_par_*scale + dm*m << " lineto\n";
	// umin
	out << "newpath\n";
	out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " moveto\n";
	out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " lineto\n";
	out << "stroke\n";

	// umax
	out << "newpath\n";
	out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " moveto\n";
	out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " lineto\n";
	out << "stroke\n";

	// vmin
	out << "newpath\n";
	out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " moveto\n";
	out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " lineto\n";
	out << "stroke\n";

	// vmax
	out << "newpath\n";
	out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " moveto\n";
	out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " lineto\n";
	out << "stroke\n";

	++mesh_it;
      }

    if(close)
      out << "%%EOF\n";
  }

  void writePostscriptMesh(const Go::Mesh2D& mesh, std::ostream &out)
  {
    writeMesh(mesh, out);
    out << "%%EOF\n";
  }

 void writePostscriptMesh(Go::LRSplineSurface& lr_spline_sf, std::ostream &out, bool highlight_gamma)
  {
    writeMesh(lr_spline_sf.mesh(), out); 
    if (highlight_gamma) highlight_non_unitary_gamma(lr_spline_sf, out);
    out << "%%EOF\n";
  }

  void writePostscriptMeshWithE(Go::LRSplineSurface& lr_spline_sf, std::ostream &out, std::unordered_set<std::pair<double,double>, LRSplineSurface::double_pair_hash>& elts_to_refine)
  {
    writeMesh(lr_spline_sf.mesh(), out); 
    out << "1 0 0 setrgbcolor\n";
    for_each(elts_to_refine.begin(), elts_to_refine.end(), [&](const std::pair<double,double>& el) {
	out << "newpath\n" << el.first << " " << el.second << " " 
	    << (lr_spline_sf.paramMax(XFIXED)-lr_spline_sf.paramMin(XFIXED))/400.0 << " " 
	    << "0 360 arc closepath\nstroke\n";
      });
    out << "%%EOF\n";
  }


  void writePostscriptMeshOverload(Go::LRSplineSurface& lr_spline_sf, std::ostream &out)
  {
    double scale = 1000;
    writeMesh(lr_spline_sf.mesh(), out, scale); 
    out << "/show-ctr {\n  dup stringwidth pop\n  -2 div 0\n  rmoveto show\n} def\n";
    for (LRSplineSurface::ElementMap::const_iterator el_it=lr_spline_sf.elementsBegin(); el_it!=lr_spline_sf.elementsEnd(); ++el_it) {
      double u = (el_it->second->umin()+el_it->second->umax())*0.5;
      double v = (el_it->second->vmin()+el_it->second->vmax())*0.5;
      int s = lr_spline_sf.basisFunctionsWithSupportAt(u, v).size();
      int r = (lr_spline_sf.degree(XFIXED)+1)*(lr_spline_sf.degree(YFIXED)+1);
      double fontsize = scale*0.5*std::min(el_it->second->umax()-el_it->second->umin(),el_it->second->vmax()-el_it->second->vmin());
      if (s > r) out << "1 0 0 setrgbcolor\n"; 
      if (s > r+1) out << "0 1 0 setrgbcolor\n"; 
      if (s > r+2) out << "0 0 1 setrgbcolor\n"; 
      if (s > r+3) out << "0 0.5 0.5 setrgbcolor\n";
      if (s > r+4) out << "0.5 0.5 0 setrgbcolor\n"; 
      if (s > r+5) out << "0.5 0 0.5 setrgbcolor\n"; 
      out << "/Times-Roman " << fontsize << " selectfont\n";
      out << scale*u << " " << scale*v-fontsize/4 << " moveto (" << s << ") show-ctr\n";
      if (s > r) out << "0 0 0 setrgbcolor\n"; 
    }
    out << "stroke\n%%EOF\n";
  }

} // end of namespace Go.
