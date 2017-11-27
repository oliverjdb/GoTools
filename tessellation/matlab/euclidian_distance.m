function [l, lx1, ly1] = euclidian_distance(P1, P2)
% 
% Euclidian distances between points in P1 and P2
%
% SYNOPSIS:
%   function [l, lx, ly] = euclidian_distance(P1, P2)
%
% DESCRIPTION:
%
% PARAMETERS:
%   P1 - Point set on form [X, Y];
%   P2 - Point set on form [X, Y];
%
% RETURNS:
%   l  - l(i,j) represents distance between Pi and Pj
%   lx1 - lx1(i,j) represents derivative of distance l_ij with respect to x1_i
%   ly1 - lx1(i,j) represents derivative of distance l_ij with respect to y1_i 
%   lx2 - (can be inferred from lx1 since it is equal to -1 * lx1)
%   ly2 - (can be inferred from ly1 since it is equal to -1 * ly1)


   %% computing distances
   [X1, X2] = ndgrid(P1(:,1), P2(:,1));
   [Y1, Y2] = ndgrid(P1(:,2), P2(:,2));
   
   dX = X1 - X2;
   dY = Y1 - Y2;
   
   l = sqrt(dX.^2 + dY.^2);
   
   %% computing derivatives
   lx1 = dX ./ l;
   ly1 = dY ./ l;