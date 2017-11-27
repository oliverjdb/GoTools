function [E, dE] = polygon_energy(poly, points, dfun, efun, R)
% @@ Assuming counterclockwise oriented polygon
% @@ assuming convex polygon, for now

   %% counting points and edges
   K = size(poly, 1);    % number of edges
      
   %% compute internal energy
   [E, dE] = interpoint_energy(points, points, dfun, efun, 'internal');
   
   %% adding boundary energy
   
   % boundary polygon energy
   [BE, dBE] = interpoint_energy(points, poly, dfun, efun, 'boundary');
   
   BOUNDARY_WEIGHT = 2; % @@ chosen because otherwise points tend to move too
                        % close to boundary points.  Perhaps this should be a
                        % tuneable parameter
   E  = E  + BOUNDARY_WEIGHT * BE; 
   dE = dE + BOUNDARY_WEIGHT * dBE;
   
   % mirror energy
   segments = [poly; poly(1,:)]; % closing polygon
   for i = 1:K
      cur_segment = segments(i:i+1,:);
  
      ind = segment_neighbor_points(cur_segment, points, R/2);

      p_ix = find(ind);
      p_ix_iy = [p_ix; p_ix + size(points, 1)];

      mpoints = mirror_points(cur_segment, points(ind,:));
      [ME, dME] = interpoint_energy(points(ind,:), mpoints, dfun, efun);
      
      E = E + ME;
      dE(p_ix_iy) = dE(p_ix_iy) + dME;
      
      % mpoints = mirror_points(cur_segment, points);
      % [ME, dME] = interpoint_energy(points, mpoints, dfun, efun);
      
      % E = E + ME;
      % dE = dE + dME;
   end
end
