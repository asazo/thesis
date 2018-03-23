function levelsets = td3d(levelsets,nt,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function levelsets = td3d(levelsets,nt,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean curvature flow of networks with arbitrary surface tensions
% in 3D. All mobilities = 1 / (surface tension).
%
% INPUT:
%   levelsets = Data structure containing initial configuration
%               of sets in the partition. Cell array, with
%               levelsets{k} = level set values for the k-th set
%               on the uniform rectangular grid (n1 x n2 array).
%               The 0-level set delineates the boundary of the
%               k-th set.
%   nt = Number of time steps to take.
%   dt = Time step size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE:
% Straight-forward, and therefore inefficient, but easy to read
% implementation of one of the algorithms from:
%
%   Esedoglu, S.; Otto, F. Threshold dynamics for networks with
%   arbitrary surface tensions. Communications on Pure and
%   Applied Mathematics. 68:5 (2015), pp. 808-864.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: Matrix of surface tensions is denoted S below. Modify
% it with the desired surface tensions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACKNOWLEDGEMENT:
%   Research leading to this software was supported by the US
%   National Science Foundation through grant DMS-0748333.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(levelsets,2); % Number of phases.
[n1,n2,n3] = size(levelsets{1});

Z = -ones(n1,n2,n3); % Aux. variable, used in calling ls2vf2D.
ind = find(Z<0);  % Aux. variable, used in calling ls2vf2D.
[sub1,sub2,sub3] = ind2sub([n1 n2 n3],ind); % Aux. variable, used in calling ls2vf2D.
indicator = 0*Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the convolution kernel KERNEL in fourier space:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = sqrt(-1);      % Imaginary number i.
wx=exp(I*2*pi/n1); % nth root of unity.
wy=exp(I*2*pi/n2); % nth root of unity.
wz=exp(I*2*pi/n3); % nth root of unity.
[x,y,z]=meshgrid([1:n1],[1:n2],[1:n3]);
x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);
KERNEL = n1*n1*(2-wx.^(x-1)-wx.^(1-x)) + n2*n2*(2-wy.^(y-1)-wy.^(1-y)) + ...
    n3*n3*(2-wz.^(z-1)-wz.^(1-z));
KERNEL = exp( -dt*KERNEL );
clear('x','y','z');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix of surface tensions:
S = ones(N,N) - eye(N,N); % Equal surface tensions.

for t=1:nt % Time stepping.

  % Comparison functions initialized:
  for k=1:N
    comparisons{k} = 0*Z;
  end

  % Convolutions:
  for k=1:N
    vals = ls2vf3D(int32(sub1),int32(sub2),int32(sub3),levelsets{k}(:),Z,n1,n2,n3);
    indicator(ind) = vals;
    convolution = real(ifftn( fftn(indicator) .* KERNEL ));
    % Convolution dropped into comparison functions it appears in:
    for ell=1:N
      if ell ~= k
        comparisons{ell} = comparisons{ell} + S(k,ell) * convolution;
      end % if.
    end % for ell.
  end % for k.
  
  % Redistribution step. Level set functions updated:
  for k=1:N
    levelsets{k} = 1e10 + 0*Z;
  end
  for k=1:N
    for ell=1:N
      if ell ~= k
        levelsets{k} = min( levelsets{k}, comparisons{ell} );
      end % if.
    end % for ell.
  end % for k.
  for k=1:N
    levelsets{k} = levelsets{k} - comparisons{k};
  end

end % for t. Time stepping ends.