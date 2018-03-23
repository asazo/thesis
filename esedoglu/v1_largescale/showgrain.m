function showgrain(grains,dims,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function showgrain(grains,dims,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = -ones(dims);
ind = grains{k,1};
[x,y,z] = ind2sub(dims,ind);
val = grains{k,2};
u(ind) = grains{k,2};
clf; isosurface(u,0.0); axis([1 dims(1) 1 dims(2) 1 dims(3)]); axis square;