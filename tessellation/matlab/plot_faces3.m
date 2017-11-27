function plot_faces3(pts, edges, faces)

   clf; hold on;
   
   % make one-based indices rather than zero-based
   edges = edges + 1;
   faces = faces + 1;
   
   for f = 1:size(faces,1)
      
      e_ixs = faces(f, 1:end-1);
      
      v = edges(e_ixs(1), :);
      
      for i = e_ixs(2:end)
         
         if v(end) == edges(i, 1)
            v = [v, edges(i, 2)];
         else
            v = [v, edges(i, 1)];
         end
      end
      
      % plot face
      v = v(1:end-1);
      if faces(f, end) < 2
         v = fliplr(v);
      end
      
      pts(v', :)
      
      patch(pts(v', 1), pts(v', 2), pts(v', 3), 'r'); view(45,45);
   end
   
   view(45, 45)
   light('position', [-1 -1 -1])
end

