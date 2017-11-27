function efun = energy_function_factory(type, radius)

   switch(type)
     case 'simple'
       efun = @(l) simple_energy_function(l, radius);
     case 'steep'
       efun = @(l) steep_energy_function(l, radius);
     otherwise 
       error('unknown energy function type.');
   end
   
end

function [e, de, de2] = simple_energy_function(l, radius)
   
   tmp = max(radius - l, 0);
   
   e = tmp.*tmp;
   e = e.*e;
   de = -4 * tmp .* tmp .* tmp;
   de2 = 12 * tmp .* tmp;
   
end

% function [e, de, de2] = simple_energy_function(l, radius)
   
%    tmp = max(radius - l, 0);
   
%    e = tmp.*tmp;
%    de = -2 * tmp;
%    de2 = 2 * (tmp > 0);
   
% end


function [e, de, de2] =  steep_energy_function(l, radius)

   epsilon = radius * 1e-4;
   r_eps = radius + eps;
   l_eps = l + eps;
   
   sqrt_e = (radius - l) ./ l_eps;
   
   e = sqrt_e .* sqrt_e;
   de = -2 * r_eps .* (radius - l).^2 ./ (r_eps.^3);
   de2 = []; % @@ unimplemented for now
    
   % setting everything outside radius to zero
   tr_ix = l>radius;   
   e(tr_ix) = 0;
   de(tr_ix) = 0;
   de2(tr_ix) = 0;
end

   