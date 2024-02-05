% function r = node_trans(D,x0,x1,rc,w)
function r = node_trans(D,rmin,rmax,rlim,rbl);
% function r = node_trans(D,varargin)

% rc = .25;
% alpha = @(D) 1./ (1 + exp(-(rc-D)./w));
% 
% r = x0 + (x1-x0) * (1 - alpha(D));


% rbl = w;
% rlim = rc;
% rmin = x0;
% rmax = x1;

% r = rmax*ones(size(D,1),1);
% r(D<rlim) = rmin;
% r((D > rc) & (D < rbl)) = rmin + (rmax-rmin) * (D((D > rc) & (D < rbl))-rlim)./(rbl-rlim);

% if D < rlim
%     r = rmin;
% elseif (D > rlim) & (D < rbl)
%     r = rmin + (rmax-rmin) * (D-rlim)./(rbl-rlim);
% else
%     r = rmax;
% end


if D < rlim
    r = rmin;
elseif (D >= rlim) && (D <= (rlim + rbl))
    r = rmin + (rmax-rmin) * (D-rlim)./rbl;
%     (D-rlim)./rbl
else
    r = rmax;
end