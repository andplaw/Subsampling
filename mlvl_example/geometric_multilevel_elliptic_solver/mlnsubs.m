% function [X,R,p] = mlnsubs(Xs,Xb,dmin,Nmin)
function [Xlvls,R,nsets] = mlnsubs(Xbg,Xb,Nmin)

Xb = node_sort(Xb,"polar");
[in,on] = inpolygon(Xbg(:,1),Xbg(:,2),Xb(:,1),Xb(:,2));
in = xor(in,on);

fc = 1.5;
K = 10;
Xtmp = Xbg;
xyb=Xb;
xyi = Xtmp(in,:);
Xlvls{1} = [xyb;xyi];
nsets=true(size(Xlvls{1},1),1);
R = [];

while size(xyb,1) >= Nmin %boundary minimum
    dir = 2;

    [xyb] = nsubs(xyb,K,fc,dir);
    [~,dbex]=knnsearch(xyb,xyb,'K',2);
    dbex = dbex(:,2);
    [idb,dii]=knnsearch(xyb,Xtmp,'K',1);
    idrem = dii <= dbex(idb)/2;

    Xtmp(idrem,:) = [];
    in(idrem) = [];
    [Xtmp,idkeep] = nsubs(Xtmp,K,fc,dir);
    in = in(idkeep);
    xyi = Xtmp(in,:);



if size(xyb,1) >= Nmin
        Xsave = [xyb;xyi];
        Xlvls{end+1} = Xsave;
        idres = ismember(Xlvls{end-1},Xsave,'rows');
        R{end+1}.id = find(idres);
        R{end}.R = ones(sum(idres),1);
        nsets(:,end+1) = ismember(Xlvls{1},Xsave,'rows');
    else
        break
    end
end

end