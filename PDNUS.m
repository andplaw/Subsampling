function xy_sub = PDNUS(xy,fc)
% Poisson Disk Inspired Non-Uniform Sample Elimination


arguments
    xy
    fc = 1.5
end

if size(xy,1)<size(xy,2)
    xy = xy';
end

rng(0);

%% algorithm

subs = zeros(size(xy,1),1);         %array to hold nodes which will remain
% xy_ = xy(randperm(length(xy)),:);   %randomize for robustness
xy_ = sortrows(xy,2);
[~,D] = knnsearch(xy_,xy_,'K',2);   %find local node sep for density param
h = fc*D(:,2)/2;                      %local density param, i.e. exclusion radius

[Idx,D] = rangesearch(xy_,xy_,max(h));  %all nodes w/in largest exclusion radius (Poisson disk)

subs(1) = 1;                        %begin with first node
j = 1;
for i = 2:length(xy_)
    idx_subs = any(subs(1:j)==Idx{i},1);    %which subsampled nodes are nearest neighbors of node i
    
    Di = D{i};
    if all(Di(idx_subs)>=(h(i)+h(idx_subs)'))||isempty(Di(idx_subs))  
        %if distance between node i and a priori subsampled nodes is more 
        %than sum of respective exlusions radii OR if no subsampled nodes
        %are even near, then
        j = j+1;        %increment how many nodes have been subsampled
        subs(j) = i;    %record node to subsample
    end
end

subs = subs(1:j);     %get rid of excess zeros
xy_sub = xy_(subs,:);   %actually subsample the nodes



end








