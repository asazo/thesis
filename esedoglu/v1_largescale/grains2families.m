function [families,famgrains] = grains2families(grains,R,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [families,famind] = grains2families(grains,R,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groups far away grains into families efficiently, for joint convolution.
% Periodic boundary conditions used.
%
% INPUT:
%   grains = Data structure (cell array) containing pixel coordinates and
%            level set function data.
%   R = Minimum distance to maintain between grains in the same family.
%   dims = Dimensions of the grid: dims = [n1 n2 n3].
% OUTPUT:
%   families = Unions of far away grains. Has same structure as grains.
%   famgrains = Which grains each family contains. Cell array, size
%               same as number of families.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global WORKSPACE LONGARRAY;

n1 = dims(1); n2 = dims(2); n3 = dims(3); % Dimension of the grid.

N = size(grains,1); % Number of grains.
famgrains = cell(100,1);

grainlabels = zeros(dims); % Grid describing which grain each pixel belongs to.

for k=1:N % Loop over grians.
  ind = grains{k,1};  % List of pixels in a nhd of grain k.
  ind2 = grains{k,2}; % Level set function val. at pixels in grain k.
  ind3 = find(ind2 > 0); % Pixels contained in the interior of grain k.
  grainlabels(ind(ind3))=k; % Populate grid of grain labels.
end

remaininggrains = zeros(1,N); % Grains not allocated to families so far.
maxlabel = 0;
families = cell(100,3); % Maximum of 100 families anticipated.

collect_ind = zeros(n1*n2*n3,1);
collect_val = zeros(n1*n2*n3,1);
collect_cval = zeros(n1*n2*n3,1);

for label=1:100 % Loop over families.

  intersects = zeros(1,N);
  collected = 0;

  listofgrains = []; % Grains contained in this family. Initialized.
  numberofgrains = 0; % Number of grains contained in this family.

  for k=1:N % Loop over the grains.
  
    if (intersects(k)==0 && remaininggrains(k)==0) % If grain k is unallocated
                                                   % and is far from other
                                                   % grains in this family
      
      % Remove grain k from remaining grains:
      remaininggrains(k) = 1;

      % Grow grain k:
      ind = grains{k,1};
      [x,y,z] = ind2sub(dims,ind);
      [x2,y2,z2] = dilation(int32(x),int32(y),int32(z),R,WORKSPACE,LONGARRAY);
      ind2 = sub2ind(dims,x2,y2,z2); % Pixels in R neighborhood of current grain.

      % Move grain to family label:
      pixingrain = size(ind,1);
      collect_ind(collected+1:collected+pixingrain) = ind;
      collect_val(collected+1:collected+pixingrain) = grains{k,2};
      collect_cval(collected+1:collected+pixingrain) = grains{k,3};
      collected = collected + pixingrain;
      
      maxlabel = label;

      % Grains ineligible for inclusion in current family due to proximity:
      gind = grainlabels(ind2); % Labels of grains found in R neighborhood.
      if max(gind) > 0
        ind2 = gind(find(gind)); % Eliminate zeros from labels.
        intersects(ind2) = 1; % Those labels marked as ineligible.
      end
      
      % Add grain k to list of grains contained in this family:
      numberofgrains = numberofgrains + 1;
      listofgrains(numberofgrains,1) = k;

    end % end if.
    
  end % (for k) Loop over grains ends.

  %families{label,1} = [families{label,1} ; ind];
  %families{label,2} = [families{label,2} ; grains{k,2}];
  families{label,1} = collect_ind(1:collected);
  families{label,2} = collect_val(1:collected);
  families{label,3} = collect_cval(1:collected);

  famgrains{label} = listofgrains;

end % (for label). Loop over families ends.

families = families(1:maxlabel,1:3);
famgrains = famgrains(1:maxlabel);