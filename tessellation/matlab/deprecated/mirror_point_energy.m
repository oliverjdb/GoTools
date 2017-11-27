function [E, dE] = mirror_point_energy(segment, points, dfun, efun)
% Return the energy associated with mirrored points across the specified segment
%
% SYNOPSIS:
%   function [E, dE] = mirror_point_energy(segment, points, efun)
%
% DESCRIPTION:
%
% PARAMETERS:
%   segment - segment across which we consider mirror points
%   points  - points that will be mirrored
%   dfun    - distance function:
%   efun    - energy function: returns the energy (and derivative)
%             associated with a distance
%
% RETURNS:
%   E  - energy value
%   dE - energy gradient
%

   mpoints = mirror_points(segment, points);
   
   %% computing distances between points and mirror points
   [l, lx, ly] = dfun(points, mpoints);
   
   %% compute distance energy E
   [e, de] = efun(l);
   
   E = sum(e(:)); % total energy
   
   %% compute partial derivatives of E
   
   dE_dx = 2 * sum(lx .* de, 2);
   dE_dy = 2 * sum(ly .* de, 2);
   
   dE = [dE_dx, dE_dy];

end
