%% Bengt's Non-Uniform Sample Elimination
clear; 
close all;

%% 
ms = 0.5; % Markersize for plots
K = 10; % Number of nearest neighbors in knnsearch
fc = 1.5; % Set factor for node elimination

ninit    = 5e+4;        % Upper bound on potential dot positions (PDPs) to use
dotmax   = 5e+5;        % Upper bound on number of dots to place
xy1 = node_drop([0 1 0 1],ninit,dotmax,@radius_trui); 

figure (1) % Show this original node set rendering
plot(xy1(:,1),xy1(:,2),'k.','MarkerSize',0.7*ms)
axis square
title(sprintf('%d xy1 dots',length(xy1)));

%%subsampling

%first subsampling

xy2 = xy1;
N = length(xy2(:,1)); % Get the number of its dots

xy2 = sortrows(xy2,2); % Sort dots from bottom and up
                       % Create nearest neighbor pointers and distances
[Idx,D] = knnsearch(xy2,xy2,'K',10);
for k = 1:N % Loop over nodes from bottom and up
    if Idx(k,1)~=0 % Check if node already eliminated
        ind = find(D(k,2:end) < fc*D(k,2));
        ind2 = Idx(k,ind+1);
        ind2(ind2<k) = []; % Mark nodes above present one, and which
        Idx(ind2,1) = 0;   % are within the factor fc of the closest one
    end                % in the original node set
end
xy3 = xy2(Idx(:,1)~=0,:); % Eliminate these marked nodes


%second subsampling
xy4 = xy3;
N = length(xy4(:,1)); % Get the number of its dots

xy4 = sortrows(xy4,2); % Sort dots from bottom and up
                       % Create nearest neighbor pointers and distances
[Idx,D] = knnsearch(xy4,xy4,'K',10);
for k = 1:N % Loop over nodes from bottom and up
    if Idx(k,1)~=0 % Check if node already eliminated
        ind = find(D(k,2:end) < fc*D(k,2));
        ind2 = Idx(k,ind+1);
        ind2(ind2<k) = []; % Mark nodes above present one, and which
        Idx(ind2,1) = 0;   % are within the factor fc of the closest one
    end                % in the original node set
end
xy5 = xy4(Idx(:,1)~=0,:); % Eliminate these marked nodes


%% plots

figure (2) % Show the subsampled image
plot(xy3(:,1),xy3(:,2),'k.','MarkerSize',1.0*ms)
axis square
title(['subsampled, ',num2str(length(xy3(:,1))),' dots']);

figure (3) % Show the twice subsampled image
plot(xy5(:,1),xy5(:,2),'k.','MarkerSize',1.0*ms)
axis square
title(['second degree subsampled, ',num2str(length(xy5(:,1))),' dots']);



figure (4) % Show both ditherings together
plot(xy1(:,1),xy1(:,2),'k.','MarkerSize',0.7*8*ms)
axis square
title(sprintf('%d xy1 dots',length(xy1))); hold on;
plot(xy3(:,1),xy3(:,2),'ko','MarkerSize',0.5*8*ms)
plot(xy5(:,1),xy5(:,2),'w.','MarkerSize',0.3*8*ms)