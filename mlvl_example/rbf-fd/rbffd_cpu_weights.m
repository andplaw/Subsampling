function [D] = rbffd_cpu_weights(Xd,X,rbfset)
%Input:
%Xd             data points (:,1:3)
%X             center/evaluation points
%rbfset         RBF settings structure containing:
%   rbfset.PHS     polyharmonic degree
%   rbfset.poly    degree of augmented polynomials
%   rbfset.op      type of operator to generate weights for ('I' interpolation or 'D' differentiation)
%   rbfset.nrat    relative stencil size (sten. size/num. poly. terms) (optional)
%   rbfset.n       stencil size (optional)

%Output:
%D.d(x,y,z)     weights for the first derivative in (x,y,z)-direction
%D.L            weights for the Laplacian operator in dim-dimensions (up to 3D)
%D.id           IDs of the neighbors that the weights should be applied to
%etc.

%Example: setup the differentiation matrix in x-direction and the Laplacian
% [x,y]=ndgrid(linspace(0,1,50),linspace(0,1,50)); %structured points (for sake of simplicity)
% Xd = [x(:) y(:)]; %data points
% X = Xd; %evaluation points (does not have to equal Xd)
% rbfset.PHS = 3; rbfset.poly = 4; rbfset.op = 'D';rbfset.nrat = 2;
% D = rbffd_cpu_weights(Xd,X,rbfset);
% dsize = size(Xd,1);
% esize = size(X,1);
% Dx=sparse(D.id1,D.id2,reshape(D.dx',[],1),esize,dsize); %Diff. matrix in x
% L=sparse(D.id1,D.id2,reshape(D.L',[],1),esize,dsize); %Laplacian

%Copyright 2022, Morten Eggert Nielsen




%Maximum number of pages for pagefun to utilize mldivide
%for systems above 64x64xNsize (will probably be fixed in new releases 2022/23)
Nsize = 682;

%Initialization of variables
D.dx=[];D.dy=[];D.dz=[];D.L=[];D.id=[];D.I=[];D.id1=[];D.id2=[];

%Number of stencils and chunks
N=size(X,1);
chunks=ceil(N/Nsize);
nset = floor(linspace(0,N,chunks+1));

%RBF-FD settings
[Ntot,dim]=size(Xd);
if isfield(rbfset,'PHS')
    m=rbfset.PHS;
    pdeg=rbfset.poly;
    if isfield(rbfset,'nrat')
        np = nchoosek(pdeg+dim,dim); %number of terms in the pth degree polynomial
        n=ceil(rbfset.nrat*np);
    elseif isfield(rbfset,'n')
        n=rbfset.n;
    else
        disp('Please define stencil size!')
    end

    Did = rbfset.op;
    if m<3
        disp(['Be aware that PHS degree is: ' num2str(m)])
    end
else
    error('Please choose a PHS-degree!')
end

%Check if the node set is large enough to hold one stencil
if Ntot<n, error('The node set is not large enough'),end

%Build KD-tree
Mdl = KDTreeSearcher(Xd);

%Loop over all chunks ("chunks of stencils")
for nchunk=1:chunks

    %Pick out the stencil centers in the current chunk
    idx=(nset(nchunk)+1):(nset(nchunk+1));
    Xe = X(idx,:);
    [Ne,~]=size(Xe);

    %NN-search
    [ID] = knnsearch(Mdl,Xe,'K',n); %NN search

    %Rearrange nearest neighbor IDs
    IDx=reshape(ID',1,n,[]);%3-dimensional array of indices of indices
    Xe = reshape((Xe'),1,[],Ne); %3-dimensional array of evaluation points


    if dim==1
        x=Xd(:,1);
        x=x(pagetranspose(IDx));
        xe = Xe(1,1,:);
        x = x-xe;
        xd = pagetranspose(x)-x;
        rd = abs(xd);
        r = abs(x);
        P = polyterms(pdeg,x); %matrix containing all polynomials up to degree "pdeg"
        LA = Lrbf(m,Did,r,(xe-xe)-x); %(xe-xe) due to centering
        LP = Lpolynml(pdeg,Did,(xe-xe)); %(xe-xe) due to centering
    elseif dim==2
        x=Xd(:,1);y=Xd(:,2);
        x=x(pagetranspose(IDx));
        y=y(pagetranspose(IDx));
        xe = Xe(1,1,:);
        ye = Xe(1,2,:);
        x = x-xe;
        y = y-ye;
        xd = pagetranspose(x)-x;
        yd = pagetranspose(y)-y;
        rd = sqrt(xd.^2 + yd.^2);
        r = sqrt(x.^2 + y.^2);
        P = polyterms(pdeg,x,y); %matrix containing all polynomials up to degree "pdeg"
        LA = Lrbf(m,Did,r,(xe-xe)-x,(ye-ye)-y); %(xe-xe) due to centering
        LP = Lpolynml(pdeg,Did,(xe-xe),(ye-ye)); %(xe-xe) due to centering
    elseif dim == 3
        x=Xd(:,1);y=Xd(:,2);z=Xd(:,3);
        x=x(pagetranspose(IDx));
        y=y(pagetranspose(IDx));
        z=z(pagetranspose(IDx));
        xe = Xe(1,1,:);
        ye = Xe(1,2,:);
        ze = Xe(1,3,:);
        x = x-xe;
        y = y-ye;
        z = z-ze;
        xd = pagetranspose(x)-x;
        yd = pagetranspose(y)-y;
        zd = pagetranspose(z)-z;
        rd = sqrt(xd.^2 + yd.^2 + zd.^2);
        r = sqrt(x.^2 + y.^2 + z.^2);
        P = polyterms(pdeg,x,y,z); %matrix containing all polynomials up to degree "pdeg"
        LA = Lrbf(m,Did,r,(xe-xe)-x,(ye-ye)-y,(ze-ze)-z); %(xe-xe) due to centering
        LP = Lpolynml(pdeg,Did,(xe-xe),(ye-ye),(ze-ze)); %(xe-xe) due to centering
    end

    %setup the full collocation matrix
    A = rd.^m;
    pterms = size(P,2); %number of terms in the augmented polynomial
    nullC=zeros(pterms,pterms,Ne); %null matrix
    C=[A P;pagetranspose(P) nullC]; %full collocation matrix

    L = [LA;LP]; %setup RHS

    %compute weights by solving Cw = L
    w=zeros(n+pterms,size(L,2),Ne);
    for ii=1:Ne
        w(:,:,ii) = C(:,:,ii)\L(:,:,ii);
    end


    if strcmp(Did,'D')
        if dim == 1
            D.dx = [D.dx;reshape(w(1:n,1,:),n,Ne,1)'];
            D.L = [D.L;reshape(w(1:n,2,:),n,Ne,1)'];
            D.id = [D.id;ID];

        elseif dim == 2
            D.dx = [D.dx;reshape(w(1:n,1,:),n,Ne,1)'];
            D.dy = [D.dy;reshape(w(1:n,2,:),n,Ne,1)'];
            D.L = [D.L;reshape(w(1:n,3,:),n,Ne,1)'];
            D.id = [D.id;ID];

        elseif dim == 3
            D.dx = [D.dx;reshape(w(1:n,1,:),n,Ne,1)'];
            D.dy = [D.dy;reshape(w(1:n,2,:),n,Ne,1)'];
            D.dz = [D.dz;reshape(w(1:n,3,:),n,Ne,1)'];
            D.L = [D.L;reshape(w(1:n,4,:),n,Ne,1)'];
            D.id = [D.id;ID];
        end
    elseif strcmp(Did,'I')
        D.I = [D.I;reshape(w(1:n,1,:),n,Ne,1)'];
        D.id = [D.id;ID];
    end

    D.id2 = [D.id2;reshape(permute(IDx, [2 1 3]),[],1)];
end

D.id1 = reshape(repmat((1:N)',1,n)',[],1);

end


function LA = Lrbf(m,op,r,varargin)
%PHS of degree m differential operators
dim=double(numel(varargin));
if strcmp(op,'D')
    if dim==1
        x=varargin{1};
        LA = [m*x.*r.^(m-2) m*m*r.^(m-2)]; %RBF RHS (Dx,DL)
    elseif dim==2
        x=varargin{1};y=varargin{2};
        LA = [x*m.*r.^(m-2) y*m.*r.^(m-2) m*m*r.^(m-2)]; %RBF RHS (Dx,Dy,DL)
    elseif dim == 3
        x=varargin{1};y=varargin{2};z=varargin{3};
        LA = [m*x.*r.^(m-2) m*y.*r.^(m-2) m*z.*r.^(m-2) m*m*r.^(m-2)]; %RBF RHS (Dx,Dy,Dz,DL)
    end

elseif strcmp(op,'I')
    LA = r.^m;
end

end

function LP = Lpolynml(d,op,varargin)
dim=double(numel(varargin)); %dimension
np = factorial(d+dim)/(factorial(d)*factorial(dim)); %number of polynomial terms

if dim==1
    x=varargin{1};
elseif dim==2
    x=varargin{1};y=varargin{2};
elseif dim == 3
    x=varargin{1};y=varargin{2};z=varargin{3};
end

N = size(x,3);

%centered stencils
if op == 'D'
    if dim==1
        LP=zeros(np,2,N); %Initialize RHSs
        if d>=1;LP(2,1,:)=1;end %Dx
        if d>=2;LP(3,2,:)=2;end %DL
    elseif dim==2
        LP=zeros(np,3,N); %Initialize RHSs
        if d>=1;LP(2,1,:)=1;LP(3,2,:)=1;end %Dx,Dy
        if d>=2;LP(5,3,:)=2;LP(6,3,:)=2;end %DL
    elseif dim==3
        LP=zeros(np,4,N); %Initialize RHSs
        if d>=1;LP(2,1,:)=1;LP(3,2,:)=1;LP(4,3,:)=1;end %Dx,Dy,Dz
        if d>=2;LP(8,4,:)=2;LP(9,4,:)=2;LP(10,4,:)=2;end %DL
    end
elseif op == 'I'
    LP=zeros(np,1,N); %Initialize RHSs
    LP(1,1,:) = 1; %I
end

end