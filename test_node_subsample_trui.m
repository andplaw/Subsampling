% --- Main script for generating the dithered trui image
close all; 
clear all;

filePath = "results/";

ms = 3;

ninit    = 5e+4;        % Upper bound on potential dot positions (PDPs) to use
dotmax   = 5e+5;        % Upper bound on number of dots to place

% --- Carry out the node dropping
xy = node_drop([0 1 0 1],ninit,dotmax,@radius_trui); 
% writematrix(xy,'trui_dithered.csv')

% --- Display the resulting dithered image
figure
plot(xy(:,1),xy(:,2),'k.','MarkerSize',ms); axis square %
% title(['Original trui image'])
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_dithered.eps','BackgroundColor','none')

A = double(imread('trui.png','PNG')); 


% figure 
% t = tiledlayout(2,2);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% nexttile % First, show the original image in regular gray scale
% image(A/255)
% title('(a) Original gray scale image')
% axis off; axis image; axis tight;
% 
% nexttile % Show the full dithered image
% plot(xy(:,1),xy(:,2),'k.','MarkerSize',markersize)
% axis square
% title('(b) Dithered version')
% % print -deps trui_01.eps

%% Subsample Yuksel

c1 = 3.44;          %coarsening factor
k = 1;          %size of density radius
alpha = 8;
r = 1;
gamma = 1.5;
beta = 0.65;

% title(t, sprintf('First subsampling, $c$=%d, $\\alpha$=%d, $k$=%.2f, $\\gamma$=%.1f, $\\beta$=%.2f',...
%                  c,alpha,k,gamma,beta),'Interpreter','latex')
[mT_y,sT_y] = timeit2(@() WNUS(xy,c1,k,alpha,r,gamma,beta),10);
xy_c_y = WNUS(xy,c1,k,alpha,r,gamma,beta); %10533 nodes
c2 = 3.1;
xy_c2_y = WNUS(xy_c_y,c2,k,alpha,r,gamma,beta); %3404 nodes


% --- Display the coarsened images alone
figure
plot(xy_c_y(:,1),xy_c_y(:,2),'k.','MarkerSize',ms); axis square %
% title(['Sub 1, WS, ', num2str(size(xy_c_y,1)), ' nodes'])
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub1_ws.eps','BackgroundColor','none')

figure
plot(xy_c2_y(:,1),xy_c2_y(:,2),'k.','MarkerSize',ms); axis square %
% title(['Sub 2, WS, ', num2str(size(xy_c2_y,1)), ' nodes'])
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub2_ws.eps','BackgroundColor','none')

% --- Display the coarsening in series
% figure
% t_y = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% plot(xy(:,1),xy(:,2),'k.','MarkerSize',ms); 
% title(sprintf('Dithered trui.png, %d nodes',...
%               size(xy,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c_y(:,1),xy_c_y(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 1, WS, %d nodes',...
%               size(xy_c_y,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c2_y(:,1),xy_c2_y(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 2, WS, %d nodes',...
%               size(xy_c2_y,1)),'Interpreter','latex')
% axis square %
% set(gcf, 'Renderer', 'painters')
% exportgraphics(t_y,filePath+'trui_subs_ws.eps','BackgroundColor','none')

%% Subsample Bengt
kb1 = 1.5101;
kb2 = 1.518;
[mT_b,sT_b] = timeit2(@() MFNUS(xy,kb1),10);
xy_c_b = MFNUS(xy,kb1);
xy_c2_b = MFNUS(xy_c_b,kb2);


% --- Display the coarsened images alone
figure % Show the subsampled image
plot(xy_c_b(:,1),xy_c_b(:,2),'k.','MarkerSize',ms)
axis square
% title(['Sub 1, MF, ',num2str(length(xy_c_b(:,1))),' nodes']);
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub1_mf.eps','BackgroundColor','none')

figure % Show the twice subsampled image
plot(xy_c2_b(:,1),xy_c2_b(:,2),'k.','MarkerSize',1.0*ms)
axis square
% title(['Sub 2, MF, ',num2str(length(xy_c2_b(:,1))),' nodes']);
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub2_mf.eps','BackgroundColor','none')

% --- Display the coarsening in series
% figure
% t_b = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% plot(xy(:,1),xy(:,2),'k.','MarkerSize',ms); 
% title(sprintf('Dithered trui.png, %d nodes',...
%               size(xy,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c_b(:,1),xy_c_b(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 1, MF, %d nodes',...
%               size(xy_c_b,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c2_b(:,1),xy_c2_b(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 2, MF, %d nodes',...
%               size(xy_c2_b,1)),'Interpreter','latex')
% axis square %
% set(gcf, 'Renderer', 'painters')
% exportgraphics(t_b,filePath+'trui_subs_mf.eps','BackgroundColor','none')

%% Subsample dexsets - bad bad not good don't use

% deg = 115;
% %fekete pts
% leja = 0;
% xy_f1 = dexsets(deg,xy,leja,[]);
% xy_f2 = dexsets(deg,xy_f1,leja,[]);
% 
% figure % Show the subsampled image
% plot(xy_f1(:,1),xy_f1(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['fekete 1, ',num2str(length(xy_f1(:,1))),' dots']);
% 
% figure % Show the twice subsampled image
% plot(xy_f2(:,1),xy_f2(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['fekete 2, ',num2str(length(xy_f2(:,1))),' dots']);
% 
% %leja pts
% leja = 1;
% xy_l1 = dexsets(deg,xy,leja,[]);
% xy_l2 = dexsets(deg,xy_l1,leja,[]);
% 
% figure % Show the subsampled image
% plot(xy_l1(:,1),xy_l1(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['leja 1, ',num2str(length(xy_l1(:,1))),' dots']);
% 
% figure % Show the twice subsampled image
% plot(xy_l2(:,1),xy_l2(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['leja 2, ',num2str(length(xy_l2(:,1))),' dots']);


%% Subsample CATCH - bad bad not good don't use

% deg = 25;
% [xy_CATCH1,~,~,~,~,~,~] = compresscub_extended(deg,xy,1./radius_trui(xy),0);
% 
% figure
% plot(xy_CATCH1(:,1),xy_CATCH1(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('First CATCH sub, %d nodes',...
%               size(xy_CATCH1,1)),'Interpreter','latex')

%% Poisson subsample for var density
kp1 = 1.4931;
kp2 = 1.5394;
[mT_p,sT_p] = timeit2(@() PDNUS(xy,kp1),10);
xy_c_p = PDNUS(xy,kp1);
xy_c2_p = PDNUS(xy_c_p,kp2);


% --- Display the coarsened images alone
figure % Show the subsampled image
plot(xy_c_p(:,1),xy_c_p(:,2),'k.','MarkerSize',1.0*ms)
axis square
% title(['Sub 1, PD, ',num2str(length(xy_c_p(:,1))),' nodes']);
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub1_pd.eps','BackgroundColor','none')

figure % Show the twice subsampled image
plot(xy_c2_p(:,1),xy_c2_p(:,2),'k.','MarkerSize',1.0*ms)
axis square
% title(['Sub 2, PD, ',num2str(length(xy_c2_p(:,1))),' nodes']);
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub2_pd.eps','BackgroundColor','none')

% --- Display the coarsening in series
% figure
% t_p = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% plot(xy(:,1),xy(:,2),'k.','MarkerSize',ms); 
% title(sprintf('Dithered trui.png, %d nodes',...
%               size(xy,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c_p(:,1),xy_c_p(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 1, PD, %d nodes',...
%               size(xy_c_p,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c2_p(:,1),xy_c2_p(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 2, PD, %d nodes',...
%               size(xy_c2_p,1)),'Interpreter','latex')
% axis square %
% set(gcf, 'Renderer', 'painters')
% exportgraphics(t_p,filePath+'trui_subs_pd.eps','BackgroundColor','none')

%% Nearest Neighbor Flagging for var density
% not tunable to any number of nodes
% subsamples are restricted to distinct numbers of nodes

% xy_c_f = NNFlag(xy);
% xy_c2_f = NNFlag(xy_c_f);
% 
% --- Display the coarsened images alone
% figure % Show the subsampled image
% plot(xy_c_f(:,1),xy_c_f(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['Sub 1, NNF, ',num2str(length(xy_c_f(:,1))),' dots']);
% 
% figure % Show the twice subsampled image
% plot(xy_c2_f(:,1),xy_c2_f(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['Sub 2, NNF, ',num2str(length(xy_c2_f(:,1))),' dots']);


%% Diversity subsampling
% from python

% fid = fopen('trui_sub1_gds.csv');
% xy_c_g = textscan(fid, '%f', 'Delimiter',',');
% xy_c_g = xy_c_g{1};
% fclose(fid);
xy_c_g = readmatrix('trui_sub1_gDS.csv');

% fid = fopen('trui_sub2_gds.csv');
% xy_c2_g = textscan(fid, '%f', 'Delimiter',',');
% xy_c2_g = xy_c2_g{1};
% fclose(fid);
xy_c2_g = readmatrix('trui_sub2_gDS.csv');


% --- Display the coarsened images alone
figure % Show the subsampled image
plot(xy_c_g(:,1),xy_c_g(:,2),'k.','MarkerSize',1.0*ms)
axis square
% title(['Sub 1, gDS, ',num2str(length(xy_c_g(:,1))),' nodes']);
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub1_gDS.eps','BackgroundColor','none')

figure % Show the twice subsampled image
plot(xy_c2_g(:,1),xy_c2_g(:,2),'k.','MarkerSize',1.0*ms)
axis square
% title(['Sub 2, gDS, ',num2str(length(xy_c2_g(:,1))),' nodes']);
ax = gca;
set(gca,'XTick',[], 'YTick', [])
exportgraphics(ax,filePath+'trui_sub2_gDS.eps','BackgroundColor','none')

% --- Display the coarsening in series
% figure
% t_g = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% plot(xy(:,1),xy(:,2),'k.','MarkerSize',ms); 
% title(sprintf('Dithered trui.png, %d nodes',...
%               size(xy,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c_g(:,1),xy_c_g(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 1, gDS, %d nodes',...
%               size(xy_c_g,1)),'Interpreter','latex')
% axis square %
% nexttile
% plot(xy_c2_g(:,1),xy_c2_g(:,2),'k.','MarkerSize',ms); axis square %
% title(sprintf('Sub 2, gDS, %d nodes',...
%               size(xy_c2_g,1)),'Interpreter','latex')
% axis square %
% set(gcf, 'Renderer', 'painters')
% exportgraphics(t_g,filePath+'trui_subs_gDS.eps','BackgroundColor','none')
% export_fig trui_subs_gDS -eps -transparent

%% Node Quality Measures

% K = 6; %average the distance of the K-1 nearest neighbors
ks = 2:14;

methods = {xy_c_y, xy_c2_y;
           xy_c_b, xy_c2_b;
           xy_c_p, xy_c2_p;
           xy_c_g, xy_c2_g};
method_names = {'Weighted';
                'Moving Front';
                'Poisson Disk';
                'generalized Diversity'};
lm = ['o';'s';'d';'x'];
ls = ["-";"--"];

CLRs = cell(size(methods,1),size(methods,2),size(ks,2),2);

for s = 1:size(methods,2)
    fig_avg = figure;
    hold on
    fig_sd = figure;
    hold on
    for m = 1:size(methods,1)
        for k = 1:size(ks,2)
            [avg,sd] = NQM(xy,methods{m,s},ks(k));
            CLRs{m,s,k,1} = avg;
            CLRs{m,s,k,2} = sd;
        end
        style = sprintf('%s%s',ls(s),lm(m));
        
        figure(fig_avg)
        p = plot(ks,[CLRs{m,s,:,1}],style,'DisplayName',method_names{m});
        set(p, 'MarkerFaceColor', get(p,'Color')); 
        title(['Average, Subsample ',num2str(s)])
        xlabel('k nearest neighbors')
        ylabel('Comparative Local Regularity')
        legend('Location','best')
        
        figure(fig_sd)
        p = plot(ks,[CLRs{m,s,:,2}],style,'DisplayName',method_names{m});
        set(p, 'MarkerFaceColor', get(p,'Color')); 
        title(['Standard Deviation, Subsample ', num2str(s)])
        xlabel('k nearest neighbors')
        ylabel('Comparative Local Regularity')
        legend('Location','best')
        
    end
    figure(fig_avg);
    ax_avg = gca();
    filename = ['CLR_avg_sub',num2str(s),'.eps'];
    exportgraphics(ax_avg,filePath+filename,'BackgroundColor','none')
    
    figure(fig_sd);
    ax_sd = gca();
    filename = ['CLR_sd_sub',num2str(s),'.eps'];
    exportgraphics(ax_sd,filePath+filename,'BackgroundColor','none')
end



%% Subsample Bengt, Rotated
fig1 = figure;
t1 = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
fig2 = figure;
t2 = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

theta=-pi/2; %radians by which to rotate
R=[cos(theta) -sin(theta); sin(theta) cos(theta)]; %Rotation matrix
Rinv=[cos(-theta) -sin(-theta); sin(-theta) cos(-theta)]; %inverse rotation matrix
for i = 0:3
    xy_rot = xy*R^(i); %initial rotate
    
    xy_rot2 = MFNUS(xy_rot,kb1); %first subsample
    xy_2 = xy_rot2*Rinv^(i); %rotate back for plot
    figure(fig1);
    nexttile
    plot(xy_2(:,1),xy_2(:,2),'k.','MarkerSize',1.0*ms)
    axis square
    title(sprintf('Rotated %d degrees',i*90))
   
    xy_rot3 = MFNUS(xy_rot2,kb2); %second subsample
    xy_3 = xy_rot3*Rinv^(i); %rotate back for plot
    figure(fig2);
    nexttile
    plot(xy_3(:,1),xy_3(:,2),'k.','MarkerSize',1.0*ms)
    axis square
    title(sprintf('Rotated %d degrees',i*90))
end
    
    

%% Swiss Cheese
%remove various circles from trui image

pts = [0.209728,0.585209;
       0.0564821, 0.947817;
       0.615107, 0.186981;
       0.72171, 0.920718;
       0.313165, 0.147505;
       0.750611, 0.617776;
       0.807,0.219];
radii = [0.06; 0.035; 0.04; 0.055; 0.035; 0.05; 0.055];

%remove circles
xy_sc = xy;
for j = 1:size(radii,1)
    idx = rangesearch(xy_sc,pts(j,:),radii(j));
    idx = idx{1};
    xy_sc(idx,:) = [];
end

figure
plot(xy_sc(:,1),xy_sc(:,2),'k.','MarkerSize',2.1); axis square %
title('Swiss Cheese trui')

figure
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

%%% Yuksel mod
c1 = 3.2286;
c2 = 3.0958;
xy_sc_c_y = WNUS(xy_sc,c1,k,alpha,r,gamma,beta);
xy_sc_c2_y = WNUS(xy_sc_c_y,c2,k,alpha,r,gamma,beta);

% --- Display the coarsened images alone
nexttile
plot(xy_sc_c_y(:,1),xy_sc_c_y(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 1, WS, %d nodes',...
              size(xy_sc_c_y,1)),'Interpreter','latex')

nexttile
plot(xy_sc_c2_y(:,1),xy_sc_c2_y(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 2, WS, %d nodes',...
              size(xy_sc_c2_y,1)),'Interpreter','latex')

% --- 


%%% Subsample Bengt
kb1 = 1.43460;
kb2 = 1.55744;
xy_sc_c_b = MFNUS(xy_sc,kb1);
xy_sc_c2_b = MFNUS(xy_sc_c_b,kb2);

nexttile
plot(xy_sc_c_b(:,1),xy_sc_c_b(:,2),'k.','MarkerSize',1.0*ms)
axis square
title(['Sub 1, MF, ',num2str(length(xy_sc_c_b(:,1))),' dots']);

nexttile
plot(xy_sc_c2_b(:,1),xy_sc_c2_b(:,2),'k.','MarkerSize',1.0*ms)
axis square
title(['Sub 2, MF, ',num2str(length(xy_sc_c2_b(:,1))),' dots']);


%%% Poisson Disk method
kp1 = 1.2884;
kp2 = 1.5409;
xy_sc_c_p = PoissonSubVarDen(xy_sc,kp1);
xy_sc_c2_p = PoissonSubVarDen(xy_sc_c_p,kp2);

% --- Display the coarsened images alone
nexttile
plot(xy_sc_c_p(:,1),xy_sc_c_p(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 1, Repulsed, PD, %d nodes',...
              size(xy_sc_c_p,1)),'Interpreter','latex')

nexttile
plot(xy_sc_c2_p(:,1),xy_sc_c2_p(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 2, Repulsed, PD, %d nodes',...
              size(xy_sc_c2_p,1)),'Interpreter','latex')

%% Repulsive Swiss Cheese

pts = [0.209728,0.585209;
       0.0564821, 0.947817;
       0.615107, 0.186981;
       0.72171, 0.920718;
       0.313165, 0.147505;
       0.750611, 0.617776;
       0.807,0.219];
radii = [0.06; 0.035; 0.04; 0.055; 0.035; 0.05; 0.055];

%remove circles
density = zeros(size(radii));
xy_sc = xy;
for j = 1:size(radii,1)
    [idx,d] = rangesearch(xy_sc,pts(j,:),radii(j));
    idx = idx{1};
    density(j) = pi*radii(j)^2/length(idx);
    xy_sc(idx,:) = [];
end

figure
plot(xy_sc(:,1),xy_sc(:,2),'k.','MarkerSize',2.1); axis square %
hold on
%plot circles
circs = cell(size(radii,1),1);
for j = 1:size(radii,1)
    th = linspace(0,2*pi,round(1e-1*2*pi*radii(j)^(2)/density(j)^(1)));
    xcirc = radii(j)*cos(th) + pts(j,1);
    ycirc = radii(j)*sin(th) + pts(j,2);
    circs{j} = horzcat(xcirc',ycirc');
    plot(xcirc,ycirc,'ro','MarkerSize',2.1)
end
%plot border
h = 0.01;
N = 1/h;
bord = vertcat(horzcat(linspace(0,1,N)',zeros(N,1)),...
               horzcat(linspace(0,1,N)',ones(N,1)),...
               horzcat(zeros(N,1),linspace(0,1,N)'),...
               horzcat(ones(N,1),linspace(0,1,N)'));
plot(bord(:,1),bord(:,2),'ro','MarkerSize',2.1)
title('Swiss Cheese trui w/ borders')

%removal
idx = rangesearch(xy_sc,bord,h/2);
xy_sc(cell2mat(cellfun(@(row) row',idx,'UniformOutput',false)),:) = [];
for j = 1:size(circs,1)
    idx = rangesearch(xy_sc,circs{j},density(j)^(1)/(2e-1*radii(j)^(1)));
    xy_sc(cell2mat(cellfun(@(row) row',idx,'UniformOutput',false)),:) = [];
end
           
%repulsion
xy_sc = repel2(xy_sc,bord,'in',10);
for j = 1:size(circs,1)
    xy_sc = repel2(xy_sc,circs{j},'out',10);
end

figure
plot(xy_sc(:,1),xy_sc(:,2),'k.','MarkerSize',2.1); axis square %
hold on
%plot circles
for j = 1:size(radii,1)
    plot(xcirc,ycirc,'ro','MarkerSize',2.1)
end
%plot border
plot(bord(:,1),bord(:,2),'ro','MarkerSize',2.1)
title('Repulsed Swiss Cheese trui')

%%% Subsampling
figure
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

%%%% Yuksel mod
xy_sc_c = WNUS(xy_sc,c1,k,alpha,r,gamma,beta);
xy_sc_c2 = WNUS(xy_sc_c,c2,k,alpha,r,gamma,beta);

% --- Display the coarsened images alone
nexttile
plot(xy_sc_c(:,1),xy_sc_c(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 1, Repulsed, WS, %d nodes',...
              size(xy_sc_c,1)),'Interpreter','latex')

nexttile
plot(xy_sc_c2(:,1),xy_sc_c2(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 2, Repulsed, WS, %d nodes',...
              size(xy_sc_c2,1)),'Interpreter','latex')

% --- 


%%%% Subsample Bengt

xy_sc2 = MFNUS(xy_sc,kb1);
xy_sc3 = MFNUS(xy_sc2,kb2);

nexttile
plot(xy_sc2(:,1),xy_sc2(:,2),'k.','MarkerSize',1.0*ms)
axis square
title(['Sub 1, Repulsed, MF, ',num2str(length(xy_sc2(:,1))),' dots']);

nexttile
plot(xy_sc3(:,1),xy_sc3(:,2),'k.','MarkerSize',1.0*ms)
axis square
title(['Sub 2, Repulsed, MF, ',num2str(length(xy_sc3(:,1))),' dots']);


%%% Poisson Disk method

xy_sc_c_p = PoissonSubVarDen(xy_sc,kp1);
xy_sc_c2_p = PoissonSubVarDen(xy_sc_c_p,kp2);

% --- Display the coarsened images alone
nexttile
plot(xy_sc_c_p(:,1),xy_sc_c_p(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 1, Repulsed, PD, %d nodes',...
              size(xy_sc_c_p,1)),'Interpreter','latex')

nexttile
plot(xy_sc_c2_p(:,1),xy_sc_c2_p(:,2),'k.','MarkerSize',ms); axis square %
title(sprintf('Sub 2, Repulsed, PD, %d nodes',...
              size(xy_sc_c2_p,1)),'Interpreter','latex')

%% Mesh Ratios

%moving front method
mr_b1 = mesh_ratio(xy_c_b);
mr_b2 = mesh_ratio(xy_c2_b);

%weighted sample elim
mr_y1 = mesh_ratio(xy_c_y);
mr_y2 = mesh_ratio(xy_c2_y);

%poisson disk
mr_p1 = mesh_ratio(xy_c_p);
mr_p2 = mesh_ratio(xy_c2_p);



%% Uniquetol() test

% xy_u = uniquetol(xy,7e-3,'ByRows',true);
% 
% figure
% plot(xy_u(:,1),xy_u(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['First sub, uniquetol(), ',num2str(length(xy_u(:,1))),' dots']);
% 
% xy_u2 = uniquetol(xy,1.45e-2,'ByRows',true);
% 
% figure
% plot(xy_u2(:,1),xy_u2(:,2),'k.','MarkerSize',1.0*ms)
% axis square
% title(['Second sub, uniquetol(), ',num2str(length(xy_u2(:,1))),' dots']);
% 
% 
% 

%% Functions 

function mr = mesh_ratio(xy)
    
    %covering radius
    %the furthest distance from a node to a vertex of its corresponding Voronoi cell
    [v,c] = voronoin(xy);
    rho = 0;
    for k = 1:size(xy,1)
        [~,D] = knnsearch(v(c{k},:),xy(k,:));
        rho = max(D(end),rho);
    end
    
    %delta
    %maxiumum distance to the nearest neighbor over whole set
    [~,D] = knnsearch(xy,xy,'K',2);
    delta = max(D(:,2)); 
    
    %mesh ratio
    mr = rho/delta; 
    
end

function [meanTime,stdTime] = timeit2(fun,N)
    
    times = zeros(N,1);
    for i = 1:N
        times(i) = timeit(fun);
    end
    meanTime = mean(times);
    stdTime = std(times);
    
end
