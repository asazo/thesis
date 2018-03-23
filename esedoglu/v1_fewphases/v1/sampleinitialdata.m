% Initial data for level sets. They can be discontinuous. 
N = 128;
f = zeros(N,N,N);
f(1:N/2,1:N/2,:) = 1;
lev0{1} = f-0.5;
f = 0*f;
f(N/2+1:N,1:N/2,:) = 1;
lev0{2} = f-0.5;
f = 0*f;
f(:,N/2+1:N,:) = 1;
lev0{3} = f-0.5;