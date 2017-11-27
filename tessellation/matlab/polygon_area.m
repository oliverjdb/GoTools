function A = polygon_area(polygon)

   A = 0;
   N = size(polygon, 1);
   origin = [0 0];
   
   for i = 1:N-1
      A = A + tri_area(origin, polygon(i,:), polygon(i+1, :));
   end
   A = A + tri_area(origin, polygon(end,:), polygon(1,:));
   
   A = abs(A);
      
end

function a = tri_area(p1, p2, p3)
   % NB: knowing that the first point is the origin, a quicker algorithm
   % would be: a = 0.5 * (-p2(2)*p3(1) + p2(1)*p3(2))
      
   a = 0.5 * det([[p1; p2; p3], ones(3,1)]);

end
