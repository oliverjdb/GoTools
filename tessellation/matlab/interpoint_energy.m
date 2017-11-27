function [E, dE] = interpoint_energy(points1, points2, dfun, efun, varargin)

   is_internal = false;
   is_boundary = false;
   if numel(varargin) > 0;
      is_internal = strcmpi(varargin{1}, 'internal');
      is_boundary = strcmpi(varargin{1}, 'boundary');
   end
      
   %% computing point distances
   [l, lx, ly] = dfun(points1, points2);
   
   if is_internal
      % remove NaNs from diagonal. They are here because point was here
      % compared against itself, so divided by zero distance (0 / 0)
      assert(all(isnan(diag(lx))));
      N = size(lx,1);
      lx(1:(N+1):N*N) = 0; % setting diagonal to zero
      ly(1:(N+1):N*N) = 0;
   end
   
   %% compute distance energy E
   [e, de] = efun(l);
   
   E = sum(e(:));
   
   %% compute partial derivatives of E
   
   fac = 2;
   if is_boundary
      % If the second point set represent the boundary polyon, these are
      % fixed, which means that derivatives should not be multiplied by two
      % (as is the case for internal point that are measured against each
      % other, as well as mirror points that are measured against mirror
      % images of each other)
      fac = 1;
   end
   
   dE_dx = fac * sum(lx .* de, 2);
   dE_dy = fac * sum(ly .* de, 2);
   
   dE = [dE_dx; dE_dy];
      
end
