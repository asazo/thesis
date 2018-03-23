function [grains] = voronoidata3d(N,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [grains] = voronoidata(N,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates initial condition with "N" grains on a grid of size "dims".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label = -1e10 * ones(dims);
n1 = dims(1); n2 = dims(2); n3 = dims(3);
W = int32(-ones(dims));
L = int32(zeros(120000000,1));

x = 1 + floor(rand(1,N)*(n1-1));
y = 1 + floor(rand(1,N)*(n2-1));
z = 1 + floor(rand(1,N)*(n3-1));

grains = cell(N,3);

% Parameter d controls the extent to which the distance funciton
% to each point in the dataset is constructed:
d = 1 + ceil( (n1*n2*n3/N)^(1/3) );
d = 2*d;

% Construct distance function, label, to the union of all the points:
for k=1:N % Loop over the random points.

  [x2,y2,z2] = meshgrid([x(k)-d:x(k)+d],[y(k)-d:y(k)+d],[z(k)-d:z(k)+d]);
  dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 + (z2-z(k)).^2 );
  x2 = 1 + mod(x2-1,n1);
  y2 = 1 + mod(y2-1,n2);
  z2 = 1 + mod(z2-1,n3);
  
  ind = sub2ind( dims , x2(:) , y2(:) , z2(:) );
  label(ind) = max(dist(:),label(ind));

end % (for k). Loop over random points ends.

% If the union of d neighborhoods of the random points do not cover
% the entire computational domain, we cannot trust the construction:
if min(label(:)) < -0.5*1e10
  error('Parameter d too small.');
end

% Associate each grid point with the random point it is closest to,
% forming the grains:
for k=1:N % Loop over the random points again.

  [x2,y2,z2] = meshgrid([x(k)-d:x(k)+d],[y(k)-d:y(k)+d],[z(k)-d:z(k)+d]);
  dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 + (z2-z(k)).^2 );
  x2 = 1 + mod(x2-1,n1);
  y2 = 1 + mod(y2-1,n2);
  z2 = 1 + mod(z2-1,n3);
  
  ind = sub2ind( dims , x2(:) , y2(:) , z2(:));
  ind2 = ind( dist(:) >= label(ind) );
  grains{k,1} = ind2;

end % (for k). Loop over random points ends.

% Dilate the grain neighborhood, and define a level set function
% on it (1 inside, -1 outside):
for k=1:N % Loop over the grains.

  ind = grains{k,1};
  [x,y,z] = ind2sub(dims,ind);
  [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),5,W,L);
  ind2 = sub2ind(dims,x2,y2,z2);
  label(ind2) = -1;
  label(ind) = 1;
  grains{k,1} = ind2;
  grains{k,2} = label(ind2);
  grains{k,3} = 0*ind2;       % Convolution vals. init to 0.

end % (for k). Loop over grains ends.
