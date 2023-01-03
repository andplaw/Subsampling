clear all
clear functions
clear
clc
close all

%% Initial Node Set
np.H = 2; %water depth
np.L = 2; %container length
np.nx = 100; %node resolution in horizontal direction
np.nz = 1.5*np.H*np.nx/np.L; %node resolution in vertical direction 
np.dx = np.L/np.nx; %node density
np.tf = 1; %simulation time
np.dt = 0.01; %time step

%Background node set params
np.xmin = 0;
np.xmax = np.L;
np.zmin = 0; %-np.H;
np.zmax = np.H; %np.H*0.5;
np.box = [-0.1*np.xmax 1.1*np.xmax -0.1*np.zmax 1.1*np.zmax];
np.ninit = 5e2;
np.dotmax = 5e4;
np.radius = @(X) norm(X)*0+np.dx;

%build background node set and boundary node set
Xold = node_drop_2d(np.box, np.ninit, np.dotmax, np.radius);

%% Coarsen?

%params
c = 3; %coarsening factor
gamma = 2; %rmin param
beta = 0.85; %rmin param
alpha = 8; %weighting param


%define current # nodes and area/vol and future # nodes
M = size(Xold,1);
N = floor(M/c);
[~,A] = boundary(Xold(:,1),Xold(:,2),0);

%set rmax
if size(Xold,2)==2
    rmax = sqrt(A/(2*sqrt(3)*N));
elseif size(Xold,2)==3
    rmax = root(A/(4*sqrt(2)*N),3);
else
    msgStr = 'Only supports 2 and 3 dimensions';
    error(msgStr)
end

%set rmin
rmin = rmax*(1-(N/M)^gamma)*beta;

%build knn model
Mdl = KDTreeSearcher(Xold);
[Idx,D] = rangesearch(Mdl,Xold,2*rmax);

%assign weights
D = cellfun(@(x) max(x,2*rmin),D,'UniformOutput',false); %Weight Limiting
W = cellfun(@(x) arrayfun(@(d) (1-d/(2*rmax))^alpha,x),D,'UniformOutput',false); %Calculate Weights

%Create Sparse Weight Matrix
I = cell2mat(cellfun(@(x) x(1)*ones(size(x,2),1),Idx,'UniformOutput',false));
J = cell2mat(cellfun(@(x) x',Idx,'UniformOutput',false));
V = cell2mat(cellfun(@(x) x',W,'UniformOutput',false));
Wij = sparse(I,J,V);

%ELIMINATE

roster = (1:M)';
counter = 1;
while counter<=M-N
    [~,j] = max(sum(Wij,2));
    Wij(j,:) = [];
    Wij(:,j) = [];
    roster(j) = [];
    counter = counter+1;
end

Xnew = Xold(roster,:);
plot(Xold(:,1),Xold(:,2),'r.',Xnew(:,1),Xnew(:,2),'gx')
figure
plot(Xold(:,1),Xold(:,2),'k.')
figure
plot(Xnew(:,1),Xnew(:,2),'k.')
        



