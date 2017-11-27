function [M, Mdx, Mdy] = distmat(points, R, varargin)

   if ~isempty(varargin) && strcmpi(varargin{1}, 'brute_force')
      [M, Mdx, Mdy] = distmat_bruteforce(points, R);
   else
      [M, Mdx, Mdy] = distmat_rapid(points, R);
   end
end

% ----------------------------------------------------------------------------

function [M, Mdx, Mdy] = distmat_bruteforce(points, R)
   
   [M, Mdx, Mdy] = euclidian_distance(points, points);
   
   too_distant_ixs = M > R;
   
   M(too_distant_ixs) = 0;
   Mdx(too_distant_ixs) = 0;
   Mdy(too_distant_ixs) = 0;
   
   M = sparse(M);
   Mdx = sparse(Mdx);
   Mdy = sparse(Mdy);
end

% ----------------------------------------------------------------------------

function [M, Mdx, Mdy] = distmat_rapid(points, R)
   
   %% Sort points into square bins of sidelength R
   [Xmin, Ymin] = deal(min(points(:,1)), min(points(:,2)));
   
   % determine x- and y index of bin each point belongs to
   points(:,1) = points(:,1) - Xmin; % translations do not affect distances
   points(:,2) = points(:,2) - Ymin; 
   
   bin_x = max(1, ceil(points(:,1)/R)); % bins are numbered from 1 and upwards
   bin_y = max(1, ceil(points(:,2)/R));
   
   % We now know that any point sorted in bin (i, j) can only have
   % sufficiently close neighbors in bins (k, l) for k=[i-1:i+1] and
   % l=[j-1:j+1].
   
   %% Assembling distance matrix
   N = size(points, 1); % total number of points
   [M, Mdx, Mdy] = deal(zeros(N));
   [bin_x_max, bin_y_max] = deal(max(bin_x), max(bin_y));
   
   for iy = 1:bin_y_max
      
      p_iy = (bin_y >= (iy-1)) & (bin_y <= (iy+1));
      
      for ix = 1:bin_x_max
         
         p_ix = (bin_x >= ix-1) & (bin_x <= ix+1);

         % measuring distances between all points in the relevant bins
         p_ixy = p_ix & p_iy; % identifies the points in the relevant bins
         p_i_glob = find(p_ixy);
         p_loc = points(p_ixy,:); % points in the relevant bins ('local' to
                                  % this iteration)
         [M_loc, M_loc_dx, M_loc_dy] = euclidian_distance(p_loc, p_loc);
         
         % Only keeping distances smaller than R
         keep = find(M_loc < R);
         
         [i_loc, j_loc] = ind2sub(size(M_loc), keep);

         ix_glob = sub2ind([N,N], p_i_glob(i_loc), p_i_glob(j_loc));
         
         M(ix_glob)   = M_loc(keep);
         Mdx(ix_glob) = M_loc_dx(keep);
         Mdy(ix_glob) = M_loc_dy(keep);
         
      end
   end

   M = sparse(M);
   Mdx = sparse(Mdx);
   Mdy = sparse(Mdy);
end
