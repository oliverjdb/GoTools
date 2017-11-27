function play_hist(polygon, x_hist, wait_time)

   clf;
   N = size(x_hist, 2);
   polygon = [polygon; polygon(1,:)];

   for i = 1:N
      clf;
      xx = reshape(x_hist(:,i), [], 2);
      plot(polygon(:,1), polygon(:,2), 'k'); hold on
      plot(xx(:,1), xx(:,2), 'b*')
      
      % tracking the first point
      plot(xx(1,1), xx(1,2), 'r*')
      pause(wait_time)
   end
   
end
