segment = [0.3 0.3; 0.7 0.7];

points = rand(80000, 2);

R = 0.2;

ind = segment_neighbor_points(segment, points, R);

figure;

plot(segment(:,1), segment(:,2), 'g-', 'linewidth', 3); hold on
plot(points(~ind, 1), points(~ind, 2), 'r.');
plot(points( ind, 1), points( ind, 2), 'b.');