function xy_sub = NNFlag(xy)

% Based on Zomolo et al paper 
% in turn based on Katz and Jameson

arguments
    xy
end

if size(xy,1)<size(xy,2)
    xy = xy';
end

rng(0);

%% 

dim = size(xy,2);
if dim==2
    nc = 6;
elseif dim==3
    nc = 13;
else
    error('Only 2 and 3 dimensions accepted')
end

%% Boundary subsampling

flags = zeros(size(xy));
k = boundary(xy(:,1),xy(:,2));

Idx_b = knnsearch(xy(k,:),xy(k,:),'K',nc);
for i = 1:length(k)
    j = k(i);
    if flags(j)==0
        flags(j) = 2;
        flags(k(Idx_b(i,2:end))) = 1;
    end
end

%% Interior subsampling

n = 1:size(xy,1);
n(k) = [];

Idx_i = knnsearch(xy(n,:),xy(n,:),'K',nc);
for i = 1:length(n)
    j = n(i);
    if flags(j)==0
        flags(j) = 2;
        flags(n(Idx_i(i,2:end))) = 1;
    end
end

%% Output

xy_sub = xy(flags==2,:);

end