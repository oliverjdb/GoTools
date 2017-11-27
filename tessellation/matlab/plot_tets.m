function plot_tets(pts, tets, color, alpha)

   hold on
   for i = 1:size(tets, 1)
      ixs = tets(i,:);
      
      ixs1 = ixs(1:3);
      ixs2 = [ixs(1:2), ixs(4)];
      ixs3 = [ixs(2:3), ixs(4)];
      ixs4 = [ixs(3), ixs(1), ixs(4)];
      
      patch(pts(ixs1,1), pts(ixs1,2), pts(ixs1,3), color, 'facealpha', alpha);
      patch(pts(ixs2,1), pts(ixs2,2), pts(ixs2,3), color, 'facealpha', alpha);
      patch(pts(ixs3,1), pts(ixs3,2), pts(ixs3,3), color, 'facealpha', alpha);
      patch(pts(ixs4,1), pts(ixs4,2), pts(ixs4,3), color, 'facealpha', alpha);
      
   end
   
end
