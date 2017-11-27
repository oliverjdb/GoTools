function [hist, corners, radius] = UnitSquare2DTest(N)
   
%P = [0.6, 0.55; 0.4 0.4; 0.9, 0.1; 0.3, 0.8];
%P = rand(N,2);
  radius = 1/N;
  tmp = linspace(1/N/2, 1-(1/N/2), N);
  [X, Y] = ndgrid(tmp ,tmp);
  P = [X(:), Y(:)];
  
  %radius = sqrt(1/N/pi);
  %radius = radius * 2;
  Nb = 2* ceil(1/radius);
  bp = [linspace(0, 1, Nb)]';
  corners = unique([repmat(0, Nb, 1), bp;
                    bp, repmat(1, Nb, 1);
                    repmat(1, Nb, 1), bp;
                    bp, repmat(0, Nb, 1)], 'rows');

  radius = 2.5 * radius;
  %corners = [ 0 0; 1 0; 1 1; 0 1];
  % corners = [0 0; 0.25 0; 0.5 0 ; 0.75 0 ; 1 0; 1 0.25; 1 0.5; 1 0.75; 1 1; ...
  %            0.75 1; 0.5 1; 0.25 1; 0 1; 0 .75; 0 .5; 0 .25];
  
  [e, d, dd] = energy(radius, P, corners);  
  hist = struct('P', P, 'err', norm(d, inf), 'e', e);

  
  f = @(u) fun(u, radius, corners);

  [v, u, history] = unitBoxBFGS(P(:), f);
    
  %function [v, u, history] = unitBoxBFGS(u0, f, varargin)
  
  
  
  
  
  
  
    
  
  
  
  
  
  small = 0.1; %0.1; %0.5;
  MAX_ITER = 5000;
  threshold = 1e-14;
  
  grad_prev = [];
  dir_prev = [];
  cur_norm = inf;
  
  for i = 1:MAX_ITER
     
     [e, d, dd] = energy(radius, P, corners);
     
     prev_norm = cur_norm;
     cur_norm = norm(d, inf);
     if cur_norm == 0;
        break
     end
     
     if prev_norm ~= cur_norm%max(abs(d)) > 1e-4
        
                
        grad_cur = -d;
        
        if isempty(grad_prev)
           grad_prev = grad_cur;
           dir_prev = grad_cur * 0;
        end
           
        gamma = ([grad_cur - grad_prev]' * grad_cur) / (grad_prev' * ...
                                                        grad_prev);
        % if (mod(i, 30) == 0)
        %   gamma = 0;
        % end
        
        dir_cur = grad_cur + gamma * dir_prev;
        
        dir_prev = dir_cur;
        grad_prev = grad_cur;
        
        %dir = -reshape(d, [], 2);
        dir = reshape(dir_cur, [], 2);

        

        
        % Identify maximum allowable step
        ix_pos = (dir)>0;
        ix_neg = (dir)<0;
        
        m1 = min(P(ix_neg) ./ abs(dir(ix_neg)));
        m2 = min( (1 - P(ix_pos)) ./ abs(dir(ix_pos)));
        if isempty(m1)
           m1 = m2;
        elseif isempty(m2)
           m2 = m1;
        end
        
        max_step = 0.9 * min(m1, m2);
        
        
        % search for internal minimum within bracket
        vals = [energy(radius, P, corners), ...
                energy(radius, P+dir*max_step, corners)];
        
        G = 0.63; %% approx golden

        t_old = max_step;
        t = max_step * (1-G);
        vt = energy(radius, P+dir*t, corners); 
        while (vt > e)
           t_old = t;
           t = t * (1-G); 
           vt = energy(radius, P+dir*t, corners);
        end
        
        
        bracket = [0, t, t_old]; % largest interval is between 
        bvals   = [e, energy(radius, P+dir*t, corners), energy(radius, P+dir * t_old, corners)];
        
        t_next = t + (t_old-t) * (1-G);
        while (bracket(3)-bracket(1)) > (1e-2 *bracket(3))
           %           fprintf('%d ; %d\n', bracket(3)-bracket(1), (1e-2 *bracket(3)));
           %           keyboard;
           vnext = energy(radius, P+dir*t_next, corners);
           if vnext < bvals(2)
              % new midpoint should be t_next
              if t_next > bracket(2)
                 bracket = [bracket(2), t_next, bracket(3)];
                 bvals = [bvals(2), vnext, bvals(3)];
                 t_next = bracket(2) + ( bracket(3) - bracket(2)) * (1-G);
              else
                 bracket = [bracket(1), t_next, bracket(2)];
                 bvals = [bvals(1), vnext, bvals(2)];
                 t_next = bracket(1) + G * (bracket(2) - bracket(1));
              end
           else
              % t_next should be new endpoint
              if t_next > bracket(2)
                 bracket = [bracket(1), bracket(2), t_next];
                 bvals = [bvals(1), bvals(2), vnext];
                 t_next = bracket(1) + G * (bracket(2) - bracket(1));
              else
                 bracket = [t_next, bracket(2), bracket(3)];
                 bvals = [vnext, bvals(2), bvals(3)];
                 t_next = bracket(2) + (1-G) * (bracket(3) - bracket(2));
              end
           end
           assert(all(diff(bracket)>0));
           
        end
        
        update = dir * bracket(3);
        
        %update = -d * small;
        
     else
        % eliminate points that do not contribute to energy
        ix = ((abs(d)<threshold) & (sum(abs(dd), 2) < threshold));
        
        if sum(ix)  == 0 
           update = -dd\d;
        else
           keep = ~ix;
           tmp = -dd(keep, keep)\d(keep);
           update = 0 * d;
           update(keep) = tmp;
        end
     end
     update = reshape(update, [], 2);
     
     P_old = P;
     P = P + update;

     ix_outside = (P<= 0) | (P >= 1);
     if sum(ix_outside(:)) > 0
        fprintf('ouch');
     end
     
     P(ix_outside) = P_old(ix_outside);
     
     % err = norm(d(~ix_outside), inf);
     err = norm(d, inf);
     hist.err = [hist.err; err];
     hist.e = [hist.e, e];

     hist.P = [hist.P; P];
     
     fprintf('Iteration %i - err: %d \n', i, err);
     
     if err < threshold
        break;
     end

     % backtrack = update;
     % count = 0;
     % while any((P(:) <= 0) | (P(:) >= 1))
     %    backtrack = backtrack / 2;
     %    P = P - backtrack;
     %    count = count + 1;
     %    if count == 1000
     %       keyboard;
     %    end
     % end
     
     
  end

  P
         
end



  
    
  function [v,g] = fun(u, radius, corners) 
     [e, d] = energy(radius, reshape(u, [], 2), corners);
     v = -e(:);
     g = -d(:);
  end
  
