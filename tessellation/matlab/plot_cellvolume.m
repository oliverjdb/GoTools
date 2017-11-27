function plot_cellvolume(mat3D, xrange, yrange, zrange, color, transparency_val)

   % entries in matrix will be treated as booleans
   hold on;
   dx = (xrange(2) - xrange(1)) / size(mat3D, 1);
   dy = (yrange(2) - yrange(1)) / size(mat3D, 2);
   dz = (zrange(2) - zrange(1)) / size(mat3D, 3);
   
   for k = 1:size(mat3D,3)
      for j = 1:size(mat3D,2)
         for i = 1:size(mat3D,1)
            if mat3D(i, j, k)
               % this box should be plotted
               hold on;
               plot_box(i, j, k, dx, dy, dz, color, transparency_val);
            end
         end
      end
   end
end

% ----------------------------------------------------------------------------
function plot_box(i, j, k, dx, dy, dz, color, alpha)
   
   %% plot left and right sides
      
   % left side
   X = repmat(dx * (i-1), 1, 4);
   Y = [dy*(j-1), dy * (j-1), dy * j, dy * j];
   Z = [dz*(k-1), dz * k, dz * k, dz * (k-1)];
   patch(X', Y', Z', color, 'facealpha', alpha);
   
   % right side
   X = repmat(dx * i, 1, 4);
   patch(X, Y, Z, color, 'facealpha', alpha);
   
   %% plot front and back sides
   
   % front side
   X = [dx * (i-1), dx * i, dx * i, dx * (i-1)];
   Y = repmat(dy * (j-1), 1, 4);
   Z = [dz*(k-1), dz*(k-1), dz*k, dz * k];
   patch(X, Y, Z, color, 'facealpha', alpha);
   
   % right side
   Y = repmat(dy * j, 1, 4);
   patch(X, Y, Z, color, 'facealpha', alpha);
   
   %% plot top and bottom sides
   X = [dx * (i-1), dx * i, dx * i, dx * (i-1)];
   Y = [dy * (j-1), dy * (j-1), dy * j, dy * j];
   Z = repmat(dz * (k-1), 1, 4);
   patch(X, Y, Z, color, 'facealpha', alpha);
   
   % right side
   Z = repmat(dz * k, 1, 4);
   patch(X, Y, Z, color, 'facealpha', alpha);
      
end
