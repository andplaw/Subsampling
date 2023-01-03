function xy = repel(xy,K,b)
%repel code

arguments
    xy
    K = 7
    b = 3
end

for c = 5:10
    [Idx,D] = knnsearch(xy,xy,'K',K);
    for n = 1:size(xy,1)
        dir = repmat(xy(n,:),K-1,1) - xy(Idx(n,2:end),:);
        dir = sum(arrayfun(@(row) row/norm(row)^b,dir),1);
        dir = dir/norm(dir);
        xy(n,:) = xy(n,:) + dir.*D(n,2)/c;
        xy(n,:) = max(vertcat(xy(n,:),[0,0]));
        xy(n,:) = min(vertcat(xy(n,:),[1,1]));
    end
end

