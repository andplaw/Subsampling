function r = node_rad(X,varargin)
persistent rfunc dmin EQ

if isempty(rfunc)
dmin = varargin{1}{1};
EQ = varargin{1}{2};
a = 1;
b=0.01*1;
% cx=0.5;cy=0.5;
cx=0.5*0;cy=0.5*0;
rfunc = @(X) a.*exp(-((cx-X(:,1)).^2+(cy-X(:,2)).^2)./b);

% a = 1;
% b=0.01*4;
% rfunc = @(X)  a.*1.0./b.^2.*exp(-((cx-X(:,1)).^2+(cy-X(:,2)).^2)./b).*(-b-cx.*X(:,1).*2.0-cy.*X(:,2).*2.0+cx.^2+cy.^2+X(:,1).^2+X(:,2).^2).*4.0;

% b=5;
% a=1/(0.5^b);
% rfunc = @(X) a * cos(b * atan2(X(:,2),X(:,1))).*(vecnorm(X,2,2).^b);


a = 1;
b=0.01*1;
% cx=0.5;cy=0.5;
cx=0.5*0;cy=0.5*0;
rdens = @(X) tanh(a*vecnorm(X-[cx cy],2,2));


% rfunc = @(X) exp(-(vecnorm(X-0.5))./b);
end
% % 1/(1+rfunc(X))
% r = dmin/(1+2*abs(rfunc(X))); % (gaussian)





% if EQ == 'laplace'
%     rbl = 0.10;
%     rc = 0.1;
%     Xc = .5;%laplace rD=0.5
%     rho2 = 3; %laplace
%     rho1 = 1; %laplace
% elseif EQ == 'poisson'
%     rbl = 0.1;
%     rc = 0.05;
%     Xc = 0;%poisson
%     rho2 = 1; %poisson
%     rho1 = 1/3; %poisson
% %         rhob = 1/4; %poisson
% 
% %     rbl = 0.3;
% %     rc = 0.1;
% %     Xc = 0;%poisson
% %     rho2 = 1; %poisson
% %     rho1 = 1/3; %poisson
% 
%     rbl = 0.10;
%     rc = 0.3;
%     Xc = 0;%laplace rD=0.5
%     rho2 = 1; %laplace
%     rho1 = 3; %laplace
% end   

% D = abs(vecnorm(X,2,2)-Xc); %(unit disk)

% 
% % %%%%%%%%% NEW
% if EQ == 'laplace'
% rbl = 0.2; %laplace
% % rbl = 0.10; %laplace
% rc = 0.2; %laplace
% % rc = 0.3; %laplace
% rho2 = 1; %laplace
% rho1 = 4; %laplace
% % rho1 = 4; %laplace
% elseif EQ == 'poisson'
% rbl = 0.20; %poisson
% rc = 0.2; %poisson
% rho2 = 1; %poisson
% rho1 = 1/4; %poisson
% end


% %%%%%%%%% PAPER
if EQ == 'laplace'
rbl = 0.20; %laplace
% rbl = 0.10; %laplace
rc = 0.3; %laplace
% rc = 0.3; %laplace
rho2 = 1; %laplace
rho1 = 3; %laplace
% rho1 = 4; %laplace
elseif EQ == 'poisson'
rbl = 0.20; %poisson
rc = 0.1; %poisson
rho2 = 1; %poisson
rho1 = 1/3; %poisson
end


% 
D=abs(vecnorm(X,2,2));
% % D = min([abs(vecnorm(X,2,2)-0) abs(vecnorm(X,2,2)-0.5)]);
% %%%%%%%%%

% D = Xc-vecnorm(X,2,2); %(unit disk)
% D = vecnorm(X-Xc,2,2);
% r = dmin*node_trans(D,rhob,rhod,rc,rc^2.0);
r = dmin*node_trans(D,rho1,rho2,rc,rbl);
end


