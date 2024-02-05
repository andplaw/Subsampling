function [L,I,R,ID,Xlvls,A0] = mlpre(Xbg,Xb,rbf_set,Nmin,plot_io)

GPU = 0;

%Node subsampling
% [Xlvls,idres,nsets] = nsubsample(Xin,Xb,dmin,Nmin);
[Xlvls,R,nsets] = mlnsubs(Xbg,Xb,Nmin);
ID.R = R; %Restriction operator by injection
ID.nsets = nsets; %Active nodes in each subset/level of nodes
ID.Xlvls = Xlvls;


% Fine low, coarse high (new)
ID.in{1} = true(size(Xlvls{1},1),1);
ID.in{1}(ismember(Xlvls{1},Xb,'rows')) = false;

rbf_set.op = 'D'; %Type of operators - 'D' for differential and 'I' for interpolation
if GPU == 1
D = rbffd_gpu_weights(Xlvls{1},Xlvls{1},rbf_set);
else
D = rbffd_cpu_weights(Xlvls{1},Xlvls{1},rbf_set);
end
% A0 = D.L;
A0 = sparse(repmat(D.id(:,1),1,size(D.L,2)),D.id,D.L,size(D.id,1),size(D.id,1));
A = D.L;
A(~ID.in{1},:) = 0;%
A(~ID.in{1},1) = 1;%Dirichlet boundaries assumed at all boundary nodes!

% A = sparse(repmat(D.id(:,1),1,size(D.L,2)),D.id,A,size(D.id,1),size(D.id,1));
% A = sparse(repmat(1:size(D.id,1),1,size(A,2)),D.id,A,size(D.id,1),size(Xlvls{1},1));
% A = sparse(D.id1,D.id2,reshape(A',[],1),size(Xlvls{1},1),size(Xlvls{1},1));



disp('Dirichlet BCs assumed on all boundaries in current implementation')
L{1}.A = A;
L{1}.id = D.id;
L{1}.id1 = D.id1;
L{1}.id2 = D.id2;

I=[];
p = size(Xlvls,2);
for i=1:p-1
    ID.in{i+1} = true(size(Xlvls{i+1},1),1);
    ID.in{i+1}(ismember(Xlvls{i+1},Xb,'rows')) = false;

    rbf_set.op = 'D'; %Type of operators - 'D' for differential and 'I' for interpolation
    if GPU == 1
    D = rbffd_gpu_weights(Xlvls{i+1},Xlvls{i+1},rbf_set);
    else
    D = rbffd_cpu_weights(Xlvls{i+1},Xlvls{i+1},rbf_set);
    end
    A = D.L;
    A(~ID.in{i+1},:) = 0;%
    A(~ID.in{i+1},1) = 1;%Dirichlet boundaries assumed at all boundary nodes!

    disp('Dirichlet BCs assumed on all boundaries in current implementation')
    L{i+1}.A = A;
    L{i+1}.id = D.id;
    L{i+1}.id1 = D.id1;
    L{i+1}.id2 = D.id2;

    rbf_setI = rbf_set;
    rbf_setI.poly = 0;
    rbf_setI.PHS = 1;
    if isfield(rbf_set,'nrat')
    rbf_setI = rmfield(rbf_setI,'nrat');
    end
    rbf_setI.n = 5;
    rbf_setI.op = 'I'; %Type of operators - 'D' for differential and 'I' for interpolation
    if GPU == 1
    Ilow = rbffd_gpu_weights(Xlvls{i+1},Xlvls{i},rbf_setI);
    else
    Ilow = rbffd_cpu_weights(Xlvls{i+1},Xlvls{i},rbf_setI);
    end
    I{i}.I = Ilow.I;
    I{i}.id = Ilow.id;



    if plot_io == 1
        if i==1,figure,end
        if size(Xlvls{i},2) == 2
            plot3(Xlvls{i}(:,1),Xlvls{i}(:,2),Xlvls{i}(:,1)*0+i,'.')
            axis equal
            hold on
            plot3(Xlvls{i}(~ID.in{i},1),Xlvls{i}(~ID.in{i},2),Xlvls{i}(~ID.in{i},1)*0+i,'o')
            plot3(Xlvls{i}(ID.in{i},1),Xlvls{i}(ID.in{i},2),Xlvls{i}(ID.in{i},1)*0+i,'s')
            xlabel('x')
            ylabel('y')
            zlabel('level')
            drawnow
        elseif size(Xlvls{i},2) == 3
            nexttile
            plot3(Xlvls{i}(:,1),Xlvls{i}(:,2),Xlvls{i}(:,3),'.')
            axis equal
            hold on
            plot3(Xlvls{i}(~ID.in{i},1),Xlvls{i}(~ID.in{i},2),Xlvls{i}(~ID.in{i},3),'o')
            xlabel('x')
            ylabel('y')
            zlabel('level')
            drawnow
        end
    end

end