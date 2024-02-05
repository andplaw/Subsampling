function [x,y,n,t] = curve2nodes(pcurve,theta_0,theta_e,dx)
%pcurve ~ either scalar input for circle or symbolic expression
if isnumeric(pcurve)
    drdphi = @(x) zeros(numel(x),1);
    pcurve = @(x) pcurve*ones(numel(x),1);
else
    drdphi = matlabFunction(diff(pcurve));
    pcurve = matlabFunction((pcurve));
end

%Initial node distribution
nr = 1e3;
theta = linspace(theta_0,theta_e,nr)'; %discretize the polar domain
rho = pcurve(theta); %discretize the polar domain
drho = drdphi(theta);
O=trapz(theta,sqrt(rho.^2 + drho.^2));

%Update uniform node distribution
nr = round(O/dx);
[px,py]=pol2cart(theta,rho);
pt = interparc(nr,px,py,'spline');
[theta,~] = cart2pol(pt(:,1),pt(:,2));
rho = pcurve(theta); %discretize the polar domain
drho = drdphi(theta);

ephi = [-sin(theta) cos(theta)]; %unit vectors in polar coordinates
er = [cos(theta) sin(theta)]; %unit vectors in polar coordinates

t = ephi.*rho + er.*drho; t = t./vecnorm(t,2,2);

n = drho.*ephi - rho.*er; n = n./vecnorm(n,2,2);
dir = -1; %outward dir=-1, inward dir=1;
t = dir*t;
n = dir*n;
[x,y] = pol2cart(theta,rho);

end

