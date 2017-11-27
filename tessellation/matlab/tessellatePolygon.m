function [points, tris, num_boundary_points] = tessellatePolygon(polygon, l)
%
% Tessellate a polygon.  Edges between elements should have approximate
% length 'l'. 
% 
% SYNOPSIS:
%   function [points, tris, boundary] = tessellatePolygon(polygon, l)
%
% DESCRIPTION:
%
% PARAMETERS:
%   polygon - polygon to tessellate.  Can be nonconvex, but should be
%             counterclockwise. 
%   l       - approximate length of internal edges (to be aimed for)
%
% RETURNS:
%   points   - the points of the final tessellation.  The boundary points are
%              listed first
%   tris     - triangles of the tessellation.  One line per triangle, and
%              each line consists of 3 indices into the points vector.
%   num_boundary_points - number of boundary points
%

%% Tessellating edges to create a unified, enclosing polygon
   
   polygon = [polygon; polygon(1,:)]; % closing polygon
   boundary_points = [];
   for i = 1:(size(polygon, 1)-1)
      cur_seg_points = tessellate_edge(polygon(i:i+1,:), l);
      boundary_points = [boundary_points; cur_seg_points(1:end-1,:)];
   end
   num_boundary_points = size(boundary_points, 1);
   
   [points, tris] = tessellate_face(boundary_points, l);
   
end
