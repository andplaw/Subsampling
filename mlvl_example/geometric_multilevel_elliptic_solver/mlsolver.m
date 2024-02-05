function [u,i,relres,usave] = mlsolver(u,f,L,I,R,ID,nu1,nu2,tol,imax)

u_ = u;

%Preallocate and initialize
i=0;
usave=zeros(imax,size(u,1));
r(1) = norm(f-sum(L{1}.A.*u(L{1}.id),2));
r0 = r(1);
rtol = tol*r0;

counter = 0;
while i < imax && r(end) > rtol % && counter<5
    i = i + 1;

    u = mlvcyc(u,f,L,I,R,ID,nu1,nu2); %ml-cycles 1 -> size(L,2)
    rnorm2 = norm(f-sum(L{1}.A.*u(L{1}.id),2));
    r(i+1) = rnorm2;
    usave(i,:) = u';
    
    if abs(log(r(i+1))-log(r(i)))<0.1
        counter = counter+1;
    else
        counter = 0;
    end

end
semilogy(r)
relres = r./r0;
usave = usave(1:i,:)';
clear mlvcyc
end



