function [xy3,idkeep] = nsubs(xy1,K,fc,dir)
% K ~ Number of nearest neighbors in knnsearch
% fc ~ Set factor for node elimination
xy2 = xy1;
N = length(xy2(:,1)); % Get the number of its dots
[xy2,I] = sortrows(xy2,dir); % Sort dots according to direction 'dir' (-1,-2,1,2)

% Create nearest neighbor pointers and distances
[Idx,D] = knnsearch(xy2,xy2,'K',K);
for k = 1:N % Loop over nodes from bottom and up
    if Idx(k,1)~=0 % Check if node already eliminated
        ind = find(D(k,2:end) < fc*D(k,2));
        ind2 = Idx(k,ind+1);
        ind2(ind2<k) = []; % Mark nodes above present one, and which
        Idx(ind2,1) = 0; % are within the factor fc of the closest one
    end % in the original node set
end
%         xy3 = xy2(Idx(:,1)~=0,:); % Eliminate these marked nodes

ID(I) = Idx(:,1);
idkeep = find(ID~=0)';
xy3=xy1(idkeep,:);

end