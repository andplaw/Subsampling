function X_bd_sorted = node_sort(X_bd,method)
    if method=="polar"
        centerx = mean(X_bd(:,1));
        centery = mean(X_bd(:,2));
        [th,~] = cart2pol(X_bd(:,1)-centerx,X_bd(:,2)-centery);
        [~,sortIdx] = sort(th,'ascend');
        
        X_bd_sorted = X_bd(sortIdx,:);

    elseif method=="knn"

        [idx,d] = knnsearch(X_bd,X_bd,'K',3);
        sortedIdx = zeros(size(idx,1),1);
        
        sortedIdx(1) = 1;
        
        for i = 1:size(idx,1)
            j = idx(i,2);
        
            if any(j == sortedIdx(:))
                k = idx(i,3);
                if any(k == sortedIdx(:))
                    if i == size(idx,1)
                        continue
                    else
                        msgStr = "Error: The algorithm got stuck";
                        error(msgStr)
                    end
                else
                    sortedIdx(i+1) = k;
                end
            else
                sortedIdx(i+1) = j;
            end
        end
        X_bd_sorted = X_bd(sortedIdx);
    end
end
