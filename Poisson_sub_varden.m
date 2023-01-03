% --- Main script for generating the dithered trui image
close all; 
clear all;

ms = 0.5;

ninit    = 5e+4;        % Upper bound on potential dot positions (PDPs) to use
dotmax   = 5e+5;        % Upper bound on number of dots to place

% --- Carry out the node dropping
xy = node_drop([0 1 0 1],ninit,dotmax,@radius_trui); 


          
          
%% Poisson Disk inspired Variable Density Subsampling
rmin = 1.15;

tic
subs = zeros(size(xy,1),1);
xy_ = xy(randperm(length(xy)),:);
[~,D] = knnsearch(xy_,xy_,'K',2);
h = rmin*D(:,2);

[Idx,D] = rangesearch(xy_,xy_,max(h));

subs(1) = 1;
j = 1;
for i = 2:length(xy_)
    idx_subs = any(subs(1:j)==Idx{i},1);
    
    Di = Di{i};
    if all(min(Di(idx_subs))>=h(i))||isempty(Di(idx_subs))
        j = j+1;
        subs(j) = i;
    end
end

subs = subs(1:j);
xy_c = xy_(subs,:);
toc

figure
plot(xy_c(:,1),xy_c(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('First sub, %d nodes',...
              size(xy_c,1)),'Interpreter','latex')







