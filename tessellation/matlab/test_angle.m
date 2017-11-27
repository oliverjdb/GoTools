p1 = [-1, 0];
p2 = [0, 0];

N = 100;
theta = linspace(0, 2*pi, N);
res = [];
dx1 = [];
dx1_alt = [];
dx2 = [];
dx2_alt = [];
dx3 = [];
dx3_alt = [];
dy1 = [];
dy1_alt = [];
dy2 = [];
dy2_alt = [];
dy3 = [];
dy3_alt = [];

for i = 1:N

   p3  = [cos(theta(i)), sin(theta(i))];
   [a, dxa, dya] = angle(p1, p2, p3);
   
   res = [res, a];
   dx1     = [dx1, dxa(1)];
   % dx1_alt = [dx1_alt, dxa_alt(1)];

   dx2     = [dx2, dxa(2)];
   % dx2_alt = [dx2_alt, dxa_alt(2)];

   dx3     = [dx3, dxa(3)];
   % dx3_alt = [dx3_alt, dxa_alt(3)];

   dy1     = [dy1, dya(1)];
   % dy1_alt = [dy1_alt, dya_alt(1)];

   dy2     = [dy2, dya(2)];
   % dy2_alt = [dy2_alt, dya_alt(2)];

   dy3     = [dy3, dya(3)];
   % dy3_alt = [dy3_alt, dya_alt(3)];
   
   
end

