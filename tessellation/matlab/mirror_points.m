function [m_points, d_m_points] = mirror_points(segment, points) 
   % 'segments' and 'points' have separate columns for x and y components.
   % 'segments' has two rows (two points); 'points' has N rows.
   % 'm_points' are the mirrored points
   % 'd_m_points' represents the transformation matrix, which can also be
   % understood as the partial derivatives of the new coordinates with
   % respect to the old ones

   %% Compute normalized direction vector and normal vector
   d = diff(segment);
   d = d / norm(d);
   
   n = [d(2), -d(1)];
      
   %% Translate system so that segment line passes through origo
   offset = segment(1,:); % first point on segment
   points = bsxfun(@minus, points, offset);
      
   %% Compute mirroring matrix (from Householder)
   M = eye(2) - 2 * (n' * n);
   
   %% Compute mirroring points
   m_points = (M * points')';
   d_m_points = M;

   %% Translate system back to original coordinates
   m_points = bsxfun(@plus, m_points, offset);
   
end
