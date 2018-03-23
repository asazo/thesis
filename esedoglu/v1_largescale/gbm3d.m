function [grains,orio] = gbm3d(nt,dt,grains,dims,ori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [grains,ori] = gbm3d(nt,dt,grains,dims,ori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D grain boundary motion with many grains. Isotropic, unequal
% surface tensions. All mobilities = 1 / (surface tension).
%
% INPUT:
%   nt = Number of time steps.
%   dt = Time step size.
%   grains = Data structure containing initial configuration of
%            grains. Cell array, with grains{k,1} = indicies of
%            pix near grain k, and grains{k,2} = level set vals
%            at those locations. Eventually, grains{k,3} contains
%            convolution values at same locations.
%   dims = Dimensions of the grid: dims=[n n] for nxn grid.
%   ori = Orientation data for grains. For now, column vector of
%         angles (2D crystallography, i.e. fiber texture).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAUTION:
%   Some C routines modify matlab data structures in place.
%   This is done in the interest of minimal memory utilization.
%   It has the undesirable side effect that input arrays cannot
%   be guaranteed to remain unaltered, especially 'grains', even
%   when the 1st output var ~= 3rd input var.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE:
%   This program is based on one of the algorithms in:
%
%   Esedoglu, S.; Otto, F. Threshold dynamics for networks with
%   arbitrary surface tensions. Communications on Pure and
%   Applied Mathematics. 68:5 (2015), pp. 808-864.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE:
%   dims = [128 128 128];
%   grains = voronoidata(1000,dims);
%   ori = 2*pi*rand(1000,1);
%   [grains,ori] = gbm3d(10,0.0005,grains,dims,ori);
%   showgrain(grains,dims,50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%   1. The order (and number) of grains listed in data structure
%      "grains" may change from one iteration to the next. But,
%      orientations in "ori" follow the grains.
%   2. Surface tensions are in subroutine updatelevelsetdata.c.
%      Currently, they are set to 1. To use e.g. Read-Shockley,
%      modify that subroutine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACKNOWLEDGEMENT:
%   Research leading to this software was supported by the US
%   National Science Foundation through grant DMS-0748333.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Auxiliary vars.
global WORKSPACE KERNEL Z LONGARRAY;

n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.
n1n2n3 = n1*n2*n3;

WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
LONGARRAY = int32(zeros(120000000,1)); % Needed for dilation.c.
Z = -ones(dims); % Another Workspace var.
                % Needed in ls2vf3D.c.

id = [1:1:size(grains,1)]'; % Used for keeping track of grains
                            % as some of them disappear.

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN TIME LOOP STARTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:nt % Main time iteration starts.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CONVOLUTION STEP:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For fast convolution, far away grains are grouped into families.
  % Then, the entire family (union of grains) is convolved at once.
  % Finally, the families are broken up into individual grains again.
  [families,famgrains] = grains2families(grains,20,dims); % Group grains into families.
  N = size(families,1); % Number of families.
  for k=1:N % Loop over families.
    cval = convolvefamily(families{k,1},families{k,2},dims); % Calulate convolution.    
    % Distribute convolution values to the grains contained in this family:
    numberofgrains = size(famgrains{k},1); % Number of grains in this family.
    listofgrains = famgrains{k}; % Column vector of grain indices.
    for ell = 1:numberofgrains % Loop over grains contained in this family.
      label = listofgrains(ell);
      ind = grains{label,1};
      grains{label,3} = cval(ind); % Read off and record the convolution vals.
    end % (for ell) Loop over grains ends.
  end % (for k) Loop over families ends.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REDISTRIBUTION STEP:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Find grains present in a nhd. of each grid pt.:
  presence = get_nhd_grains(grains,dims(1)*dims(2)*dims(3));
  % Redistribution according to Esedoglu-Otto paper:
  updatelevelsetdata(presence,grains,id,ori);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REMOVE EMPTY GRAINS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Not necessary, but speeds up the algorithm.
  [grains,id] = removeemptygrains(grains,dims,id);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFRESH GRAIN BUFFERS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Level set data for each grain extends beyond the interface,
  % to a tubular neighborhood. That neighborhood has to be
  % updated once the interface has (potentially) moved.
  N = size(grains,1); % Number of grains.
  for k=1:N % Loop over grains.
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    cval = grains{k,3}; % Convolution vals. at those pixels.
    Z(ind) = val;      % Lev. set. representation on grid.
    posind = ind(val>0); % Pixels in the interior of grain.
    [x,y,z] = ind2sub(dims,posind);
    [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),5,WORKSPACE,LONGARRAY); % Dilation.
    ind2 = sub2ind(dims,x2,y2,z2);
    val2 = Z(ind2);    % Level set vals.
    Z(ind2) = -1;      % Reset Z back to all -1's.
    Z(ind) = cval - 1; % Convolution values - 1.
    cval2 = Z(ind2);   % Convolution vals - 1.
    Z(ind2) = -1;      % Reset Z back to all -1's.
    grains{k,1} = ind2;   % Refresh grain's data structure.
    grains{k,2} = val2;   % Ditto.
    grains{k,3} = cval2 + 1; % Ditto.
  end % (for k). Loop over grains ends.

end % (for t). Main time iteration ends.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN TIME LOOP ENDS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Order of grains may have changed; return new orientation array.
N = size(grains,1); % Number of grains.
for k=1:N % Loop over grains.
  orio(k) = ori(id(k));
end

end