function plot_edges3(pts, edges)

%   clf; hold on;
   
   for e = (edges+1)'
      plot3(pts(e, 1), pts(e, 2), pts(e,3), '-gx');
   end   
   view(45, 45);
end

