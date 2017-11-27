function plot_tri(pts, tri, color, alpha)

   hold on
   
   for i = 1:size(tri, 1)
      
      ixs = [tri(i,:), tri(i,1)]';
      
      if size(pts,2) == 2
         plot(pts(ixs,1), pts(ixs,2), ['*-' color]);
      else
         %plot3(pts(ixs,1), pts(ixs,2), pts(ixs, 3), '*-r');
         plot3(pts(ixs,1), pts(ixs,2), pts(ixs, 3), ['-' color]);
         patch(pts(ixs,1), pts(ixs,2), pts(ixs,3), color, 'facealpha', alpha, ...
               'edgecolor', 'k', 'edgealpha', alpha);
      end
   end
   
end

