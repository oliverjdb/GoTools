function [x, iter, x_hist] = cgrad2D(fun, x0, polygon)
%
% Prototype conjugate gradient minimization on bounded 2D domain
%
% SYNOPSIS:
%   function x = cgrad2D(fun, x0, polygon)
%
% DESCRIPTION:
%
% PARAMETERS:
%   fun     - takes x and returns value and gradient
%   x0      - initial point positions [xcoords; ycoords]
%   polygon - bound in 2D (assumed convex for now)
%
% RETURNS:
%   x - optimized points

   MAX_ITER = 200; % 2000
   threshold = 1e-5;%sqrt(eps);
   
   prev_grad = [];
   prev_dir = [];
   x = x0(:);
   N = numel(x)/2; % @ for large number of points, perhaps should be set to a
                   % smaller number
   x_hist = x;
   
   for i = 1:MAX_ITER
      
      [f, df] = fun(x);
      
      cur_norm = norm(df, Inf);
      if cur_norm < threshold
         break;
      end

      %% Identify new search diraction
      cur_grad = -df;
      
      if isempty(prev_grad)
         prev_grad = cur_grad;
         prev_dir = cur_grad * 0;
      end
      
      gamma = max(((cur_grad - prev_grad)' * cur_grad) / ...
                  (prev_grad' * prev_grad), 0);

      if mod(i, N) == 0
         gamma = 0; % reset
      end
      
      cur_dir = cur_grad + gamma * prev_dir;
      prev_dir = cur_dir; % should not be normalized
      prev_grad = cur_grad;
      cur_dir = cur_dir / norm(cur_dir); % for numerical stability below
      
      %% find max steplength
      fac = 0.5; %0.95; % step slightly less than max to avoid crashing on the boundary
      max_step = fac * identify_max_step(x, cur_dir, polygon);
      
      %% Search for internal minimum within bracket
      step = linmin(fun, x, cur_dir, max_step, [f, fun(x + cur_dir * max_step)]);
      
      x = x + step * cur_dir;
      x_hist = [x_hist, x];
   end
   iter = i;
   
end

% ----------------------------------------------------------------------------

function mstep = identify_max_step(x, stepdir, polygon)
% @@ implemented very inefficiently for now
   x = reshape(x(:), [], 2);
   stepdir = reshape(stepdir(:), [], 2);
   N = size(x,1);
   pointwise_max_step = zeros(N, 1);
   
   for i = 1:N
      pointwise_max_step(i) = max_step_for_point(x(i,:), stepdir(i,:), polygon);
   end
   
   mstep = min(pointwise_max_step);
   
end

% ----------------------------------------------------------------------------

function mstep = max_step_for_point(p, stepdir, polygon)

   K = size(polygon, 1);
   closed_polygon = [polygon; polygon(1,:)];
   tol = sqrt(eps);
   sdn = norm(stepdir);
   if (sdn < eps)
      mstep = Inf;
      return;
   end
   
   stepdir_norm = stepdir / sdn;
   
   for i = 1:K
      % p + s * stepdir = t * A + (1-t) * B
      [~, t, s] = ray_intersect(p, stepdir_norm, closed_polygon(i:i+1, :));
      if (s >=0) && (t >= -tol) && (t<=1+tol)
         mstep = s / sdn;
         return;
      end
   end
   error('we should never get here');
end
% ----------------------------------------------------------------------------

function [ip, t, s] = ray_intersect(p, stepdir, seg)

   A = seg(1,:)';
   B = seg(2,:)';
   
   rhs = p(:) - B;
   M = [A-B, -stepdir(:)];
   
   if abs(det(M)) < sqrt(eps)
      % we consider zero determinant, and return infinite values
      t = Inf;
      s = t;
      ip = [];
      return;
   end
   
   tmp = M\rhs;
   t = tmp(1);
   s = tmp(2);
   ip = t * A + (1-t) * B;
end

% ----------------------------------------------------------------------------

function step = linmin(fun, x, stepdir, max_step, fvals)
% @@ line minimization by golden search.  brute force method
   
   G = 0.63; % golden ratio
   
   %% Find bracket
   t_old = max_step;
   t = max_step * (1-G);
   f_old = fvals(2);
   ft = fun(x + t * stepdir);
   while (ft > fvals(1))
      t_old = t;
      f_old = ft;
      t = t * (1-G);
      ft = fun(x + t * stepdir);
   end
   
   bracket = [0, t, t_old];
   bvals = [fvals(1), ft, f_old];
   
   t_next = t + (t_old - t) * G;
   
   while(bracket(3) - bracket(1)) > (1e-2 * bracket(3))
      f_next = fun(x + stepdir * t_next);
      if f_next < bvals(2)
         % new midpoint should be t_next
         if t_next > bracket(2)
            bracket = [bracket(2), t_next, bracket(3)];
            bvals = [bvals(2), f_next, bvals(3)];
            t_next = bracket(2) + ( bracket(3) - bracket(2)) * (1-G);
         else
            bracket = [bracket(1), t_next, bracket(2)];
            bvals = [bvals(1), f_next, bvals(2)];
            t_next = bracket(1) + G * (bracket(2) - bracket(1));
         end
      else
         % t_next should be new endpoint
         if t_next > bracket(2)
            bracket = [bracket(1), bracket(2), t_next];
            bvals = [bvals(1), bvals(2), f_next];
            t_next = bracket(1) + G * (bracket(2) - bracket(1));
         else
            bracket = [t_next, bracket(2), bracket(3)];
            bvals = [f_next, bvals(2), bvals(3)];
            t_next = bracket(2) + (1-G) * (bracket(3) - bracket(2));
         end
      end
      assert(all(diff(bracket)>0));
   end      
   step = bracket(3);
end
