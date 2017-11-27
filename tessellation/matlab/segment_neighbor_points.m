function ind = segment_neighbor_points(segment, points, R)
% identify points within distance R from the (inside) side of the segment, and
% not outside the two lines perpendicular to the segment and going through each
% of the endpoints.  Returns a vector of zeros and ones, identifying which
% points are fulfilling those criteria

%% Computing normalized tangent vector
   t = diff(segment);
   t = t/norm(t);
   
%% Cross product criterion (ensure that points are within distance R)
   
vec1 = bsxfun(@minus, points, segment(1,:));  
vec2 = bsxfun(@minus, points, segment(2,:));  
      
cprod = -vec1(:,1) * t(2) + vec1(:,2) * t(1); % vec1 or vec2; doesn't matter here
      
%% Scalar product criterion

sprod1 = vec1(:,1) * t(1) + vec1(:,2) * t(2);
sprod2 = vec2(:,1) * t(1) + vec2(:,2) * t(2);

ind = (cprod > 0) & (cprod < R) & (sprod1 > 0) & (sprod2 < 0);

%ind = logical(ones(size(ind))); % @@
end