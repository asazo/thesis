function showlevelsets(lev)

N = size(lev,2);
if N<1
  error('No level sets.');
end

if ndims(lev{1}) == 2
  for k=1:N
    contour(lev{k},[0 0],'b');
    hold on
  end
  hold off
  axis square
end

if ndims(lev{1}) == 3
  [n1,n2,n3] = size(lev{1});
  for k=1:N
    isosurface(lev{k},0);
    hold on
  end
  hold off
  axis([1 n1 1 n2 1 n3]);
end