function cval = convolvefamily(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cval = convolvefamily(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolves the characteristic function of union of grains whose
% level set representation is stored in "ind" and "val" input vars
% with the kernel "KERNEL" (global variable, assumed present).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global KERNEL Z; % These global variables *must* have been already
                 % allocated before calling this function! Both are
                 % of dimension "dims".

%% Convert level set data to volume of fluid representation:
[x,y,z] = ind2sub(dims,ind);
vf = ls2vf3D(int32(x),int32(y),int32(z),val,Z,dims(1),dims(2),dims(3));

%% Carry out the convolution:
Ku = zeros(dims);
Ku(ind) = vf; % Characteristic function of the union of grains.
Ku = real(ifftn( fftn(Ku) .* KERNEL ));
%cval = Ku(ind); % Convolution values, returned.
cval = Ku;