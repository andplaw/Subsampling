function [u,v] = mlcyc(k,u,f,L,I,R,ID,nu1,nu2)
persistent L1

if isempty(L1)
    L1 = sparse(L{1}.id1,L{1}.id2,reshape(L{1}.A',[],1),size(L{1}.id,1),size(L{1}.id,1));
end
gamma = 1;


% % The recursive linear multigrid algorithm cf. Trottenberg (multigrid book)
% %Pre-smoothing
% [u] = relax_scheme(u,f,L{k}.A,L{k}.id,nu1);
% 
% %Coarse grid correction
% % d = f-sum(L{k}.A.*u(L{k}.id),2); %Compute defect
% % d = d(R{k}.id); %Restrict defect
% d = f(R{k}.id,1)-sum(L{k}.A(R{k}.id,:).*u(L{k}.id(R{k}.id,:)),2); %Compute and restrict defect
% 
% if k==2
%     v = L1\d; %Compute approximate solution by direct solver
% elseif k > 2
%     for gg=1:gamma
%         v = mlcyc(k-1,d*0,d,L,I,R,ID,nu1,nu2);
%     end
% end
% v = sum(I{k-1}.I.*v(I{k-1}.id),2); %Interpolate the correction
% u = u + v; %Compute the corrected approximation on k-grid
% % u = u + sum(I{k-1}.I.*v(I{k-1}.id),2);
% 
% %Post-smoothing
% [u] = relax_scheme(u,f,L{k}.A,L{k}.id,nu2);
% 
% end


% %% The recursive linear multigrid algorithm (Saad book)

if k == 1
       u = L1\f; %Compute approximate solution by direct solver
else
    
%Pre-smoothing
[u] = relax_scheme(u,f,L{k}.A,L{k}.id,nu1);

%Coarse grid correction
% d = f-sum(L{k}.A.*u(L{k}.id),2); %Compute defect
% d = d(R{k-1}.id); %Restrict defect
d = f(R{k}.id,1)-sum(L{k}.A(R{k}.id,:).*u(L{k}.id(R{k}.id,:)),2); %Compute and restrict defect


for gg=1:gamma
    v = mlcyc(k-1,d*0,d,L,I,R,ID,nu1,nu2);
end
% v = sum(I{k}.I.*v(I{k}.id),2); %Interpolate the correction
% u = u + v; %Compute the corrected approximation on k-grid
u = u + sum(I{k-1}.I.*v(I{k-1}.id),2);

%Post-smoothing
[u] = relax_scheme(u,f,L{k}.A,L{k}.id,nu2);

end