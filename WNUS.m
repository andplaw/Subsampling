function Xcoarse = WNUS(Xfine,c,k,alpha,r,gamma,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input: array, n-by-2 or n-by-3
%
%  Output: array, n-by-2 or n-by-3
%
%--------------------------------------------------------------------------
%!
%! \original author: Cem Yuksel
%! \derivative author: Andrew Lawrence
%! 
%!
%! \brief  Implementation of the weighted sample elimination method.
%!
%! This file includes an implementation of the weighted sample elimination 
%! method for generating Poisson disk sample sets.
%!
%! Blue noise (Poisson disk) sample sets produce high-quality sampling. They
%! often lead to lower noise and better convergence with Monte Carlo sampling.
%! They provide a uniform sample distribution over a sampling domain. Unlike
%! regular random sampling, Poisson disk sample sets avoid placing any two
%! samples too close together (determined by a Poisson disk radius).
%! 
%! The weighted sample elimination method implemented in this file generates a 
%! subset of samples with blue noise (Poisson disk) characteristics from a given 
%! input sample set. The weighted sample elimination method is simple, 
%! computationally efficient, and suitable for any sampling domain. It produces 
%! high-quality blue noise sample sets with a relatively large average Poisson 
%! disk radius without the need for specifying a Poisson disk radius. It also 
%! allows progressive (adaptive) sampling and it is efficient for high-
%! dimensional sampling. However, it does not guarantee maximal coverage.
%!
%! More details can be found in the original publication:
%!
%! Cem Yuksel. 2015. Sample Elimination for Generating Poisson Disk Sample Sets. 
%! Computer Graphics Forum 34, 2 (May 2015), 25-32. 
%! http://www.cemyuksel.com/research/sampleelimination/
%!
%-------------------------------------------------------------------------------
%
% Copyright (c) 2016, Cem Yuksel <cem@cemyuksel.com>
% All rights reserved.
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy 
% of this software and associated documentation files (the "Software"), to deal 
% in the Software without restriction, including without limitation the rights 
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
% copies of the Software, and to permit persons to whom the Software is 
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all 
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
% SOFTWARE.
% 
%--------------------------------------------------------------------------

% params
arguments
    Xfine           %Initial node set
    c = 4;          %coarsening factor
    k = 2;          %size of density radius
    alpha = 8;      %weighting param
    r = 1;          %rmax coeff
    gamma = 1.5;    %rmin param
    beta = 0.65;    %rmin param
end

%define current # nodes and area/vol and future # nodes
M = size(Xfine,1);
N = floor(M/c);

%build knn model
Mdl = KDTreeSearcher(Xfine);

%find max distance between nearest neighbors
[~,d] = knnsearch(Mdl,Xfine,'K',2);
dmax = max(d(:,2));
clear d

%find all points within k*dmax
[Idx,D] = rangesearch(Mdl,Xfine,k*dmax);
icnt = cell2mat(cellfun(@(x) size(x,2),Idx,'UniformOutput',false)); %count nearest neighbors

%set rmax based on local point density
if size(Xfine,2)==2
    A = pi*(k*dmax)^2;
    rmax = r*sqrt(c*A./(2*sqrt(3).*icnt));
elseif size(Xfine,2)==3
    A = 4/3*pi*(k*dmax)^3;
    rmax = r*root(c*A./(4*sqrt(2).*icnt),3);
else
    msgStr = 'Only supports 2 and 3 dimensions';
    error(msgStr)
end
clear icnt
% rmax = (max(rmax)+min(rmax))-rmax;

%set rmin
% rmin = rmax.*(1-(N/M)^gamma)*beta;

% rmax = (rmax+rmin)/2+c/3*rmax; %+0.75*rmax
% rmin = rmax.*(1-(N/M)^gamma)*beta;
rmax = rmax.*((1+beta-beta*(1/c)^gamma)/2+c/3);

% figure
% X = Xfine(:,1);
% Y = Xfine(:,2);
% Z = rmax;

% F = TriScatteredInterp(X, Y, Z);
% %the steps of tx and ty will allow you change smoothness
% tx = min(X):0.01:max(X);
% ty = min(Y):0.01:max(Y);
% [qx,qy] = meshgrid(tx,ty);
% qz = F(qx,qy);
% mesh(qx,qy,qz);
% hold on
% Z = rmin;
% F = TriScatteredInterp(X, Y, Z);
% %the steps of tx and ty will allow you change smoothness
% qz = F(qx,qy);
% mesh(qx,qy,qz);

rmax = num2cell(rmax);
% rmin = num2cell(rmin);

%manually do a rangesearch with rmax array
rmax_cnt = cellfun(@(x,r) sum(x<=r),D,rmax,'UniformOutput',false); %count points in radius rmax
Idx = cellfun(@(x,cnt) x(1:cnt),Idx,rmax_cnt,'UniformOutput',false); %remove points outside of rmax radius
D = cellfun(@(x,cnt) x(1:cnt),D,rmax_cnt,'UniformOutput',false); %remove points outside of rmax radius
clear rmax_cnt

%assign weights
% D = cellfun(@(x,rmin_) max(x,2*rmin_),D,rmin,'UniformOutput',false); %Weight Limiting
% clear rmin
W = cellfun(@(x,rmax_) arrayfun(@(d) (1-d/(2*rmax_))^alpha,x),D,rmax,'UniformOutput',false); %Calculate Weights
clear rmax D

%Create Sparse Weight Matrix
I = cell2mat(cellfun(@(x) x(1)*ones(size(x,2),1),Idx,'UniformOutput',false));
J = cell2mat(cellfun(@(x) x',Idx,'UniformOutput',false));
clear Idx
V = cell2mat(cellfun(@(x) x',W,'UniformOutput',false));
clear W
Wij = sparse(I,J,V);
clear I J V

% figure
% X = Xfine(:,1);
% Y = Xfine(:,2);
% Z = full(sum(Wij,2));
% 
% F = TriScatteredInterp(X, Y, Z);
% %the steps of tx and ty will allow you change smoothness
% tx = min(X):0.01:max(X);
% ty = min(Y):0.01:max(Y);
% [qx,qy] = meshgrid(tx,ty);
% qz = F(qx,qy);
% mesh(qx,qy,qz);

%ELIMINATE
roster = (1:M)';
% moo(M) = struct('cdata',[],'colormap',[]);
% figure
% i = 0;
while size(roster,1)>N
%     i = i+1;
    [~,j] = max(sum(Wij,2));
    Wij(j,:) = [];
    Wij(:,j) = [];
    roster(j) = [];

%     plot(Xfine(roster,1),Xfine(roster,2),'k.','MarkerSize',0.5); axis square %
%     drawnow
%     moo(i) = getframe;
end

Xcoarse = Xfine(roster,:);

end