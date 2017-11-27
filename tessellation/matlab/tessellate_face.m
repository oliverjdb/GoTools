function [points, tris] = tessellate_face(boundary_points, l)
%
% Subdivide the planar face enclosed by the polygon specified by
% 'boundary_points' into a set of triangles.
% 
% SYNOPSIS:
%   function [points, tris] = tessellate_face(boundary_points, l)
%
% DESCRIPTION:
%
% PARAMETERS:
%   boundary_points - Boundary points specifying the enclosing polygon.  The
%                     polygon may be nonconvex, but the points should be
%                     presented in counterclockwise order
%   l               - Target distance between each point in the
%                     triangulation.  Actual distance may deviate from this.
%
% RETURNS:
%   points - List of points, starting with the boundary points
%   tris   - list of triangles, each triangle is specified as three indices
%            into the points vector
%
   
   ipoints = determine_interior_points(boundary_points, l);
   
   tris = []%;triangulate_points(boundary_points, ipoints);
   
   points = [boundary_points; ipoints];
   
   
end

% ----------------------------------------------------------------------------

function ipoints = determine_interior_points(bpoints, l)
   
   points = init_startpoints(bpoints, l);

   R = 2 * l; % will ensure a reasonable amount of points in the neighborhood
                % of each point
   
   efun = energy_function_factory('simple', R);
   dfun = @euclidian_distance;
   
   [ipoints, iter, hist] = ...
       cgrad2D(@(x) polygon_energy(bpoints, reshape(x(:), [], 2), dfun, efun, R), ...
               points(:), bpoints);
   
   %@@
   iter 
   
   ipoints = reshape(ipoints, [], 2);
end

% ----------------------------------------------------------------------------

function spoints = init_startpoints(bpoints, l)

   poly_area = polygon_area(bpoints);

   % estimate the number of points that would cover the area with an average
   % interpoint distance of 'l':
   N = 3 * poly_area / (pi * l^2);
   
   % computing bounding box
   lx = max(bpoints(:,1)) - min(bpoints(:,1));
   ly = max(bpoints(:,2)) - min(bpoints(:,2));
   
   box_area = lx * ly;
   
   box = [min(bpoints(:,1)), min(bpoints(:,2)); ...
          max(bpoints(:,1)), max(bpoints(:,2))];

   % If we want N points to fall within the polygon, how much must we fill up
   % the box area? 
   N_box = ceil(N * box_area / poly_area);
   
   nx = ceil(sqrt(N_box * lx / ly));
   ny = ceil(N_box/nx);
   
   dx = lx / (nx+1);
   dy = ly / (ny+1);
   
   xpts = linspace(box(1,1) + dx/2, box(2,1) - dx/2, nx);
   ypts = linspace(box(1,2) + dy/2, box(2,2) - dy/2, ny);
   
   [px, py] = ndgrid(xpts, ypts);
   box_points = [px(:), py(:)]; 
   
   [in, on] = inpolygon(box_points(:,1), box_points(:,2), ...
                        bpoints(:,1), bpoints(:,2));
   
   spoints = box_points(in & ~on, :);

   % perturb points a little bit, to help kick off the optimization
   perturb = rand(size(spoints)) * min(dx, dy) * 2e-1; %@@ useful?
   spoints = spoints + perturb;

end

% ----------------------------------------------------------------------------

function triangulate_points(bpoints, ipoints)
   
   
   
end
