function [u] = mlvcyc(u,f,L,I,R,ID,nu1,nu2)
persistent e r Lp

p = size(L,2);
if isempty(e)
    e = cell(1,p); 
    for j = 2:p
        e{j} = zeros(size(R{j-1}.id,1),1);
    end
    r=e;

    Lp = sparse(repmat(L{p}.id(:,1),1,size(L{p}.A,2)),L{p}.id,L{p}.A,size(L{p}.id,1),size(L{p}.id,1));
end


%Pre-smoothing
u = relax_scheme(u,f,L{1}.A,L{1}.id,nu1);
% u = relax_schemeB(u,f,B,L{1}.A,L{1}.id,nu1);

%Compute and restrict defect/residual
r{2} = f(R{1}.id,1)-sum(L{1}.A(R{1}.id,:).*u(L{1}.id(R{1}.id,:)),2);%Compute defect/residual
% r{2}(B(R{1}.id)) = f(B(R{1}.id),1)-sum(L{1}.A(B(R{1}.id),:).*u(L{1}.id(B(R{1}.id),:)),2);%Compute defect/residual

for j=2:p-1
    %Pre-smoothing
    e{j} = relax_scheme(r{j}*0,r{j},L{j}.A,L{j}.id,nu1);
    % e{j} = relax_schemeB(r{j}*0,r{j},B(R{j-1}.id),L{j}.A,L{j}.id,nu1);

    %Compute and restrict defect/residual
    r{j} = r{j}-sum(L{j}.A.*e{j}(L{j}.id),2); %Compute defect/residual
    r{j+1} = r{j}(R{j}.id); %Restrict defect/residual

    % r{j}(B(R{j-1}.id)) = r{j}(B(R{j-1}.id))-sum(L{j}.A(B(R{j-1}.id),:).*e{j}(L{j}.id(B(R{j-1}.id),:)),2); %Compute defect/residual
    % if any(r{j}(~B(R{j-1}.id)))
    %     error('residual values have leaked outside of the interior subdomain')
    % end
    % r{j+1}(B(R{j}.id)) = r{j}(B(R{j}.id)); %Restrict defect/residual

end
%Compute coarse grid solution by direct solver
e{p} = Lp\r{p};
% e{p}(B(R{p-1}.id)) = Lp(B(R{p-1}.id),B(R{p-1}.id))\(r{p}(B(R{p-1}.id))-Lp(B(R{p-1}.id),~B(R{p-1}.id))*e{p}(~B(R{p-1}.id)));


for j=p-1:-1:2
    %Interpolate the correction
    e{j} = e{j} + sum(I{j}.I.*e{j+1}(I{j}.id),2);
    % intermediate = sum(I{j}.I.*e{j+1}(I{j}.id),2);
    % e{j}(B(R{j-1}.id)) = e{j}(B(R{j-1}.id)) + intermediate(B(R{j-1}.id));

    %Post-smoothing
    e{j} = relax_scheme(e{j},r{j},L{j}.A,L{j}.id,nu2);
    % e{j} = relax_schemeB(e{j},r{j},B(R{j-1}.id),L{j}.A,L{j}.id,nu2);
end

%Interpolate the correction
u = u + sum(I{1}.I.*e{2}(I{1}.id),2);
% intermediate = sum(I{1}.I.*e{2}(I{1}.id),2);
% u(B) = u(B) + intermediate(B);

%Post-smoothing
u = relax_scheme(u,f,L{1}.A,L{1}.id,nu2);
% u = relax_schemeB(u,f,B,L{1}.A,L{1}.id,nu2);

end