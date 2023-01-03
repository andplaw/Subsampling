function [avg,sod] = NQM(xy0,xy1,K)
    
    arguments
        xy0
        xy1
        K = 2
    end
    
    Asub = all(ismembertol(xy0,xy1,1e-10),2); %subsampling index
    
    %--- nearest neighbor distance(average) metric
    h0 = NND(xy0,K); %for initial node set
    h0 = h0(Asub);
    h1 = NND(xy1,K); %for comparison node set
    
    %--- norm of error
    avg = norm((h0-h1));
    
    if K>2
        %--- nearest neighbor distance (standard of deviation) metric
        g0 = NNSD(xy0,K); %for initial node set
        g0 = g0(Asub);
        g1 = NNSD(xy1,K); %for comparison node set

        %--- norm of error
        sod = norm((g0-g1));
    else
        sod = NaN;
    end
    
    
    function h = NND(xy,K)
        [~,D] = knnsearch(xy,xy,'K',K);
    
        if K>2
            h = mean(D(:,2:end),2);
        else
            h = D(:,2);
        end
        
        %Normalize distribution
        h = (h-min(min(h)))/max(max(h));
    end
    
    function g = NNSD(xy,K)
        [~,D] = knnsearch(xy,xy,'K',K);
    
        if K>2
            g = std(D(:,2:end),0,2);
            %Normalize distribution
            g = (g-min(min(g)))/max(max(g));
        else
            msgStr = 'The standard deviation of one point does not exist'; 
            error(msgStr)
        end
    end
    
end