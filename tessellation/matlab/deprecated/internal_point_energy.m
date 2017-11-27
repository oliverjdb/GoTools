function [E, dE] = internal_point_energy(points, dfun, efun)

   %% computing distances between points
   [l, lx, ly] = dfun(points, points);
   
   % fix the diagonal elements, which should all be NaNs here
   assert(all(isnan(diag(lx))));
   assert(all(isnan(diag(ly))));
   lx(isnan(lx)) = 0; % @@ this is assuming there are no other NaNs
   ly(isnan(ly)) = 0; % @@ this is assuming there are no other NaNs
   
   %% computing distance energy E
   [e, de] = efun(l);
   
   E = sum(e(:));
   
   %% compute partial derivatives of E

   dE_dx = 2 * sum(lx .* de, 2);
   dE_dy = 2 * sum(ly .* de, 2);
   
   dE = [dE_dx, dE_dy];
   
end
