function [val, der, der2] = energy(radius, points, boundary_points)

   % computing inter-point energies
   [ival, ider, ider2] = interpoint_energy(radius, points, boundary_points);
   
   %ival = 0; ider = 0; ider2 = 0;
   % adding boundary energies
   [bval, bder, bder2] = boundary_energy(radius, points);
   
   val  = ival  + 1 * bval;
   der  = ider  + 1 * bder;
   der2 = ider2 + 1 * bder2;
   
end

function [val, der, der2] = boundary_energy(radius, points)

   [ex1, dex1, ddex1] = dist_energy_2D(2 * points(:,1)      , radius);
   [ex2, dex2, ddex2] = dist_energy_2D(2 * (1 - points(:,1)), radius);
   [ey1, dey1, ddey1] = dist_energy_2D(2 * points(:,2)      , radius);
   [ey2, dey2, ddey2] = dist_energy_2D(2 * (1 - points(:,2)), radius);
   
   val = sum(ex1 + ex2 + ey1 + ey2);
   
   der = 2 * [dex1 - dex2; dey1 - dey2];
   
   der2 = diag(4 * [ddex1 + ddex2; ddey1 + ddey2]);
   
end
   
   
function [val, der, der2] = interpoint_energy(radius, points, boundary_points) 

   N = size(points, 1); % number of interior points
   all_points = [points; boundary_points]; 
   
   [X, Y] = ndgrid(all_points(:,1), all_points(:,2));
   Y = Y'; 
   
   dX = X - X'; % x-differences between points
   dY = Y - Y'; % y-differences between points

   D = sqrt(dX.^2 + dY.^2); % entry (i,j) : distance between point p_i and p_j
   
   [evals, devals, ddevals] = dist_energy_2D(D, radius);
   
   val = 0.5 * sum(evals(:));
   
   Dm = D + eye(size(D,1)); % modified D to avoid 0/0 in expressions below
   Edx = devals .* dX ./ Dm; 
   Edy = devals .* dY ./ Dm;
   
   der = [sum(Edx, 2), sum(Edy, 2)];
   der = der(1:N, :);
   der = der(:);

   dxx = -1./(Dm.^2) .* (ddevals .* (dX.^2) + devals .* (D - (dX.^2)./Dm));
   dxx(logical(eye(size(dxx,1)))) = -sum(dxx, 2);
   dxx = dxx(1:N, 1:N);
   
   dyy = -1./(Dm.^2) .* (ddevals .* (dY.^2) + devals .* (D - (dY.^2)./Dm));
   dyy(logical(eye(size(dyy, 1)))) = -sum(dyy, 2);
   dyy = dyy(1:N, 1:N);
   
   dxy = dX .* dY ./ (Dm.^2) .* (devals./Dm - ddevals);
   dxy(logical(eye(size(dxy, 1)))) = -sum(dxy, 2);
   dxy = dxy(1:N, 1:N);

   der2 = [dxx, dxy; dxy, dyy];
           
end

function [val, der, der2] = dist_energy_2D(dist, radius)

   %@@ Test: using a two times as long radius for boundary points
   % radius = dist*0+radius; % matrix of size 'dist'
   % radius(N+1:end, :) = 2 * radius(1,1);
   % radius(:, N+1:end) = 2 * radius(1,1);
   
   r = radius - dist;
   
   if size(r,1) == size(r,2)
      % square matrix : we are use the function to compute point distances
      r(logical(eye(size(r)))) = 0; % no energy associated with a point and
                                    % itself
   else
      % do nothing.  We are computing distance to boundary
   end
   
   val  = r.*r;
   der  = -2 * r;
   der2 = 0*r+2;
   
   val (r<=0) = 0;
   der (r<=0) = 0;
   der2(r<=0) = 0;
   
end

