function [u] = relax_scheme(u,f,L,ID,numit)
% Relaxation solver script
% u ~ unknowns/initial guess (N x 1)
% f ~ rhs (N x 1)
% L ~ (N x n) weight matrix cf. discretization method
% ID ~ global ID's for the 'n' neighbors to each of the 'N' equations
% numit ~ number of iterations

id = ID(:,2:end);
T = L(:,2:end)./L(:,1);
c = f./L(:,1);

%Gauss-Seidel
for j=1:numit
    for i=1:size(id,1)
%                     u(i,1) = c(i) - T(i,:)*u(id(i,:),1);
        u(i,1) = c(i,1);
        for k=1:size(id,2)
            u(i,1) = u(i,1) - T(i,k)*u(id(i,k));
        end
    end
end