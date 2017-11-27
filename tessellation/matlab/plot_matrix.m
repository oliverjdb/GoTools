function plot_matrix(mat, xrange, yrange)

   c = 'brgy';
   mmin = min(mat(:));
   mat = mat - mmin + 1; % matrix has values from 1 upwards
   hold on;
   dx = (xrange(2) - xrange(1)) / size(mat, 1);
   dy = (yrange(2) - yrange(1)) / size(mat, 2);
   for i = 0:size(mat,1)-1
      for j = 0:size(mat,2)-1
         X = [xrange(1) + i * dx, ...
              xrange(1) + (i+1) * dx, ...
              xrange(1) + (i+1) * dx, ...
              xrange(1) + i * dx, ...
              xrange(1) + i * dx];
              
         Y = [yrange(1) + j * dy, ...
              yrange(1) + j * dy, ...
              yrange(1) + (j+1) * dy, ...
              yrange(1) + (j+1) * dy,...
              yrange(1) + j * dy];
              
         patch(X, Y, c(mat(i+1, j+1)), 'edgecolor', 'k');
      end
   end
   
end


