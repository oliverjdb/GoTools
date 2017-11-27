function pts = tessellate_edge(segment, l)
%
% Tessellate (subdivide) the edge into smaller parts, of approximately length
% 'l'.  If the segment is shorter than l, it will be returned unchanged.  If
% it is longer, at least two segments will be returned.
% 
% SYNOPSIS:
%   function pts = tessellate_edge(segment, R)
%
% DESCRIPTION:
%
% PARAMETERS:
%   segment - line segment, specified by a start point and end point
%   l       - the aimed-for distance between each point on the segment
%
% RETURNS:
%   pts - the points that specify the subdivided segment.  The first and last
%   point will be the original endpoints of the segment.
%

   L = norm(segment(1,:) - segment(2,:));
   
   % total number of segments to be returned
   num_segs = ceil(L/l);

   % actual length of each returned segment
   l_seg = L/num_segs;
   
   % number of interior points
   num_int_points = num_segs - 1;
   
   % computing position of interior points
   t = cumsum(repmat(l_seg/L, num_segs, 1));
   t = t(1:end-1); 
   
   if numel(t) > 0
      ipoints  = bsxfun(@times, segment(1,:), (1-t)) + ...
                 bsxfun(@times, segment(2,:), t);
   else
      ipoints = [];
   end
   
   % assembling result and returnin
   pts = [segment(1,:); ipoints; segment(2,:)];
         
end
