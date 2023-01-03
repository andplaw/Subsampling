function xy1 = repel2(xy1,xy2,orientation,K,b)
%repel code

arguments
    xy1     %node set to repel
    xy2     %node set of boundary
    orientation %are xy1 nodes inside of or outside of xy2 boundary nodes
    K = 7
    b = 3
end

if strcmp(orientation,'in')
    o = true;
elseif strcmp(orientation,'out')
    o = false;
end

for c = [5,7,10]
    [Idx,D] = knnsearch(xy1,xy2,'K',K);
    for n = 1:size(xy2,1)
        xy3 = xy1(Idx(n,:),:);
        for m = 1:size(xy3,1)
            dir = repmat(xy3(m,:),K-1,1) - xy3([1:m-1,m+1:size(xy3,1)],:);
            dir = sum(arrayfun(@(row) row/norm(row)^b,dir),1);
            dir = dir/norm(dir);
            xy1(Idx(n,m),:) = xy1(Idx(n,m),:) + dir.*D(n,2)/c;
            if o~=inpolygon(xy1(Idx(n,m),1),xy1(Idx(n,m),2),xy2(:,1),xy2(:,2))
                xy1(Idx(n,m),:) = xy1(Idx(n,m),:) - dir.*D(n,2)/c;
            end
%             xy1(Idx(n,m),:) = max(vertcat(xy1(n,:),[0,0]));
%             xy1(Idx(n,m),:) = min(vertcat(xy1(n,:),[1,1]));
        end
    end
end
