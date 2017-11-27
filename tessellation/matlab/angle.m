function [gamma, ddx, ddy] = angle(p1, p2, p3)
   
   % computing angle
   d1 = p2(:) - p1(:);
   d2 = p3(:) - p2(:);
   
   d1_norm = sqrt(d1' * d1);
   d2_norm = sqrt(d2' * d2);
   
   dcross = d1(1) * d2(2) - d1(2) * d2(1);
   
   gamma = sign(dcross) * acos(d1'*d2 / (d1_norm * d2_norm));
   
   % computing gradient
   
   dcross_inv = 1/dcross;
   dprod_inv = 1/(d1'*d2);
   
   % the region around gamma=0 should use a asin-based formula to avoid singularity
   if abs(gamma) < pi/4
      [ddx, ddy] = sin_based_derivative(dprod_inv, d1, d2, d1_norm, d2_norm, gamma);
   else
      [ddx, ddy] = cos_based_derivative(dcross_inv, d1, d2, d1_norm, d2_norm, gamma);
   end
   
end

function [ddx, ddy] = cos_based_derivative(dcross_inv, d1, d2, d1_norm, d2_norm, gamma)
   
   ddx(1) = dcross_inv * ( d2(1) - d1'*d2 / (d1_norm^2) * d1(1));
   ddx(3) = dcross_inv * ( -d1(1) + d1'*d2 / (d2_norm^2) * d2(1));
   ddx(2) = -(ddx(1) + ddx(3));
   
   ddy(1) = dcross_inv * (  d2(2) - d1'*d2 / (d1_norm^2) * d1(2));
   ddy(3) = dcross_inv * ( -d1(2) + d1'*d2 / (d2_norm^2) * d2(2));
   ddy(2) = -(ddy(1) + ddy(3));

   % ddx(1) = dcross_inv * ( (p3(1) - p2(1)) + d1'*d2 / (d1_norm^2) * (p1(1) - p2(1)));
   % ddx(3) = dcross_inv * ( (p1(1) - p2(1)) + d1'*d2 / (d2_norm^2) * (p3(1) - p2(1)));
   % ddx(2) = -(ddx(1) + ddx(3));
   
   % ddy(1) = dcross_inv * ( (p3(2) - p2(2)) + d1'*d2 / (d1_norm^2) * (p1(2) - p2(2)));
   % ddy(3) = dcross_inv * ( (p1(2) - p2(2)) + d1'*d2 / (d2_norm^2) * (p3(2) - p2(2)));
   % ddy(2) = -(ddy(1) + ddy(3));
   
end

function [ddx, ddy] = sin_based_derivative(dprod_inv, d1, d2, d1_norm, d2_norm, gamma)

   % computing gradient, alternative approach (based on asin rather than acos)
   
   ddx(1) = dprod_inv * ( -d2(2) + sin(gamma) * d2_norm/d1_norm * d1(1) );
   ddx(3) = dprod_inv * ( d1(2) - sin(gamma) * d1_norm/d2_norm * d2(1) );
   ddx(2) = -(ddx(1) + ddx(3));
   
   ddy(1) = -dprod_inv * ( -d2(1) + sin(gamma) * d2_norm/d1_norm * d1(2) );
   ddy(3) = dprod_inv * ( d1(1) - sin(gamma) * d1_norm/d2_norm * d2(2) );
   ddy(2) = -(ddy(1) + ddy(3));

   % ddx_alt(1) = dprod_inv * ( (p2(2) - p3(2)) + sin(gamma) * d2_norm/d1_norm * (p2(1) - p1(1)) );
   % ddx_alt(3) = dprod_inv * ( (p2(2) - p1(2)) + sin(gamma) * d1_norm/d2_norm * (p2(1) - p3(1)) );
   % ddx_alt(2) = -(ddx_alt(1) + ddx_alt(3));
   
   % ddy_alt(1) = -dprod_inv * ( (p2(1) - p3(1)) + sin(gamma) * d2_norm/d1_norm * (p2(2) - p1(2)) );
   % ddy_alt(3) = dprod_inv * ( (p2(1) - p1(1)) + sin(gamma) * d1_norm/d2_norm * (p2(2) - p3(2)) );
   % ddy_alt(2) = -(ddy_alt(1) + ddy_alt(3));

end
