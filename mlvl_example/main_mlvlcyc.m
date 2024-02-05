%multilevel solver test script
warning('off')
addpath('geometric_multilevel_elliptic_solver')
addpath(genpath('misc scripts'))
addpath('rbf-fd')
clc
clear
close all
SAVE = 0;
PLOT = 0;

for kkk=1 %1:2
    close all
    if kkk==1
        EQ = 'laplace';
        Nmin = 1*120;
    elseif kkk==2
        EQ = 'poisson';
        Nmin = 60;
    end
    
    for iii=1 %1:4 %paper polynomials
        % for iii = 1

        if EQ == 'laplace'
            Ntotal = 1e3*69*2.^(0:4); %1e3*70*2.^(0:4);
            % Ntotal = 1e3*70*2.^(0:2);
        elseif EQ == 'poisson'
            Ntotal = 1e3*70*2.^(0:4);
            % Ntotal = 1e3*70*2.^(0:2);
        end
        
        for jjj=1:numel(Ntotal)
            clear functions
            [iii jjj]
            % Test Linear Multi-Grid solver lmg for solving Poisson's equation on unit square
            dim = 2;
            N = Ntotal(jjj);
            imax = 100;

            %%%Structured boundary nodes (unit square)
            nxy = round(sqrt(N));
%             [x,y]=ndgrid(linspace(0,1,nxy));
%             Xb = [x(:) y(:)];
%             rho_b = x(2)-x(1);
            rho_b = 1/sqrt(Ntotal(jjj)) * 1.0;
            if EQ == 'poisson'
                rho_b = rho_b*3*1;
            end


            IB=1
            if IB==1
                syms t
                %Generate node sets (background, dirichlet, neumann and interior)
                dscale = 1;
                %     rD = 0.25 * (2 + sin(2*t) - 0.01*cos(5*t-pi/2) + 0.63*sin(6*t-0.1)); %butterfly
                rD = 0.5; %unit disk
                [xx,yy,nn,~] = curve2nodes(rD,0,2*pi,rho_b*(1/dscale));
                xD=xx(1:end-1,1);yD=yy(1:end-1,1);nn=nn(1:end-1,:);
                Xb=[xD(:) yD(:)]*dscale + 0.0;
            end



            %Scatter nodes and repel
            ninit = nxy; % Upper bound on potential dot positions (PDPs) to use
            dotmax = 5e+8; % Upper bound on number of dots to place
            radii = @node_rad;
            % radii = @(X,varargin) rho_b;
            Xs = node_drop([-0.1 1.1 -0.1 1.1]-0.5, 5e4, 5e6,radii,rho_b,EQ);
            [~,Xs] = repel(Xs,Xb,Xb,@(xx,yy) rho_b); %with or without repelling
            Xbg = Xs;
            
            %RBF-FD configuration
            rbf_set.poly = 2*iii; %polynomial degree
            rbf_set.PHS = 3;
            rbf_set.nrat = 2;


            %% Multilevel pre-computation (operators, identifiers and node sets)
            tic
            [L,I,R,ID,Xlvls,A0] = mlpre(Xbg,Xb,rbf_set,Nmin,PLOT);
            Xin = Xlvls{1}(ID.in{1},:);
            Xd = Xlvls{1};
            toc


            %% Setup test problem
            if EQ == 'laplace'
                % % Laplace - oscillatory
                b=10;
                a=1/(0.5^b);
                usol = @(a,b,X) a * cos(b * atan2(X(:,2),X(:,1))).*(vecnorm(X,2,2).^b); %Strang page 281
                rhs = @(a,b,X) zeros(size(X,1),1);
            elseif EQ == 'poisson'
                % % %Poisson - gaussian
                a = 2;
                b=0.01*1;
                cx=0.5*0;cy=0.5*0;
                rhs=@(a,b,X)a.*1.0./b.^2.*exp(-((cx-X(:,1)).^2+(cy-X(:,2)).^2)./b).*(-b-cx.*X(:,1).*2.0-cy.*X(:,2).*2.0+cx.^2+cy.^2+X(:,1).^2+X(:,2).^2).*4.0;
                usol = @(a,b,X)a.*exp(-((cx-X(:,1)).^2+(cy-X(:,2)).^2)./b);
            end
            
            %% Initialize source and solution vectors
            f = [usol(a,b,Xb); rhs(a,b,Xin)];
            u0 = 0*[usol(a,b,Xb); usol(a,b,Xin).*rand(size(Xin,1),1)*0];
            uorig = usol(a,b,Xlvls{1});
              
            %% Geometric multilevel solver
            tol = 1e-16;
            tic
            [u,i,relres,usave] = mlsolver(u0,f,L,I,R,ID,2,1,tol,imax);
            cpu_ml(jjj)=toc
            
            
            %% Store data
            rel_res{jjj} = relres
            q{jjj} = relres(2:end)./relres(1:end-1);
            maxrelerr{jjj} = max(abs(usave-uorig))./max(abs(uorig)); %absolute
            abserr(jjj) = maxrelerr{jjj}(end);
            max(abs(usave(:,1)-uorig))./max(abs(uorig))

            %% Compute solution using direct solver
            if 1 == 1
                tic
                A = sparse(L{1}.id1,L{1}.id2,reshape(L{1}.A',[],1),size(L{1}.A,1),size(L{1}.A,1));
                uu=A\f;
                cpu_direct(jjj)=toc
                cpu_speedup(jjj) = cpu_direct(jjj)/cpu_ml(jjj)
                [max(abs(u-uorig)) max(abs(uu-uorig))]
            else
                max(abs(u-uorig))
            end

            dsave(jjj) = rho_b;
            Nsave(jjj) = size(Xlvls{1},1)
            Xsave{jjj} = Xlvls{1};
            Esave{jjj} = abs(u-uorig)./max(abs(uorig));
            Nlvls{jjj} = size(Xlvls,2);
        end

        %% PLOT FUNCTIONALITIES


        if iii==1
            if 1==1
                % Solution and error plots
                figure
                ms = 3;
                % figure
                scatter3(Xsave{end}(:,1)-0.5*0,Xsave{end}(:,2)-0.5*0,usol(a,b,Xsave{end}),ms,usol(a,b,Xsave{end}),'filled'),
                % scatter3(Xsave{end}(:,1)-0.5*0,Xsave{end}(:,2)-0.5*0,u,ms,u,'filled'),
                % axis equal
                % cbh = colorbar ; %Create Colorbar
                % cbh.Ticks = linspace(0, 1, 5);
                % cbh.Limits = [0 1];
                % cbh.FontSize = 16;
                ax = gca;
                % ax.XAxis.FontSize = cbh.FontSize;
                % ax.YAxis.FontSize = ax.XAxis.FontSize;
                % ax.ZAxis.FontSize = ax.XAxis.FontSize;
                ax.LineWidth = 2;
                xlim([-0.52 0.5])
                xticks([-0.5  0  0.5])
                ylim([-0.52 0.5])
                yticks([-0.5  0  0.5])
                % zlim([0 1])
                if EQ == 'poisson'
                    %     zticks([0 0.25 0.50 0.75 1])
                    zticks([0:0.5:2])
                elseif EQ == 'laplace'
                    zticks([-1:0.5:1])
                end
                ax.XAxis.FontSize = 16;
                ax.YAxis.FontSize = 16;
                ax.ZAxis.FontSize = 16;
                %  colorbar
                grid off
                xlabel('$x$','interpreter','latex','fontsize',25)
                ylabel('$y$','interpreter','latex','fontsize',25)
                zlabel('$u(x,y)$','interpreter','latex','fontsize',25)
                set(gcf,'Position',[680 608.2000 766.6000 489.8000])
                if SAVE == 1
                    savefig(gcf,['ml0_' EQ '.fig'])
                    saveas(gca, ['ml0_' EQ '.eps'],'epsc');
                end
            end

            if 1==1
                figure
                tiledlayout(1,size(Xlvls,2),'TileSpacing','Compact');
                % t = tiledlayout(1,size(Xlvls,2),'TileSpacing','Compact');
                % ms = 0.0;
                for ii=1:size(Xlvls,2)
                    ms = 0.5+ii
                    %     ms = ms+0.5
                    nexttile
                    plot(Xlvls{ii}(:,1),Xlvls{ii}(:,2),'.k','MarkerSize',ms), grid off
                    axis equal off tight
                    xlabel('x','interpreter','latex'),
                    ylabel('y','interpreter','latex'),
                    xlim([0 0.55]),ylim([0 0.55])
                end
                set(gcf,'Position', [1 539.4000 1916 579.2000])
                % fontsize(gcf, 14, 'points')
                if SAVE == 1
                    savefig(gcf,['ml1_' EQ '.fig'])
                    saveas(gca, ['ml1_' EQ '.eps'],'epsc');
                end
            end

        end




        %% Performance plots

        all_marks = {'o','s','d','^','v','>','<','p','h','*','.','x','+'};
        % all_colors = {"r","g",'b','c','m'};
        blue = [0 0.4470 0.7410];
        red = [0.8500 0.3250 0.0980];
        yellow = [0.9290 0.6940 0.1250];
        purple = [0.4940 0.1840 0.5560];
        green = [0.4660 0.6740 0.1880];
        cyan = [0.3010 0.7450 0.9330];
        all_colors = {blue, red, cyan, purple, green, yellow};

        Nlegend = [];
        ms = 6;
        addsze = 2;
        ticksize = 16+addsze;
        lblsize = 20+addsze;
        lgndsize = 18;
        Lwidth = 1.5;
        pLwidth = 1.5;
        subpos = [1:4;5:8;9:12;13:16];

        figure
        for jjj=1:size(rel_res,2)
            % figure(200)
            % figure
            % subplot(4,4,subpos(iii,1))
            subplot(1,4,1),
            semilogy(1:size(rel_res{jjj},2),rel_res{jjj},'linewidth',pLwidth,'color',all_colors{jjj},'marker',all_marks{jjj},'markersize',ms,'MarkerFaceColor',all_colors{jjj});hold on;grid on
            ax = gca;
            ax.XAxis.FontSize = ticksize;
            ax.YAxis.FontSize = ticksize;
            xlim([0 imax])
            xticks([0:10:imax])
            set(gca,'LineWidth',Lwidth)
            ylim([1e-16 1e1])
            yticks([1e-16 1e-12 1e-8 1e-4 1e-0])
            xlabel('\textbf{Iterations}','interpreter','latex','fontsize',lblsize)
            ylabel('\textbf{Relative residual}','interpreter','latex','fontsize',lblsize)
            % legend([num2str(Nsave(:))],'Location','best')

            % subplot(4,4,subpos(iii,2))
            subplot(1,4,2)
            semilogy(1:size(maxrelerr{jjj},2),maxrelerr{jjj},'linewidth',pLwidth,'color',all_colors{jjj},'marker',all_marks{jjj},'markersize',ms,'MarkerFaceColor',all_colors{jjj}),hold on,grid on
            ax = gca;
            ax.XAxis.FontSize = ticksize;
            ax.YAxis.FontSize = ticksize;
            xlim([0 imax])
            xticks([0:10:imax])
            set(gca,'LineWidth',Lwidth)
            xlabel('\textbf{Iterations}','interpreter','latex','fontsize',lblsize)
            ylabel('\textbf{Maximum relative error}','interpreter','latex','fontsize',lblsize)
            ylim([1e-13 1e1])
            yticks([1e-12 1e-8 1e-4 1e-0])

            % subplot(4,4,subpos(iii,3))
            subplot(1,4,3),
            Nfunc = 2*abserr(1)*((1./sqrt(Nsave)).^(rbf_set.poly))./((1./sqrt(Nsave(1))).^(rbf_set.poly));
            % Nfunc = abserr(1)*(dsave.^(rbf_set.poly))./(dsave(1).^(rbf_set.poly));
            loglog(1./sqrt(Nsave),abserr,'ok-',1./sqrt(Nsave),Nfunc,'--k','linewidth',pLwidth,'markersize',ms,'MarkerFaceColor','k'),grid on
            set(gca,'XDir','reverse')
            ax = gca;
            ax.XAxis.FontSize = ticksize;
            ax.YAxis.FontSize = ticksize;
            set(gca,'LineWidth',Lwidth)
            xlabel('\mbox{\boldmath $\rho_{mean}= 1/\sqrt{N}$} ','interpreter','latex','fontsize',lblsize)
            ylabel('\textbf{Maximum relative error}','interpreter','latex','fontsize',lblsize)
            xlim([1e-3 1e-2])
            xticklabels({'0.001','0.01'})
            ylim([1e-13 1e1])
            yticks([1e-12 1e-8 1e-4 1e-0])

            % subplot(4,4,subpos(iii,4))
            subplot(1,4,4)
            loglog(Nsave,cpu_ml./cpu_ml(1),'ko-',Nsave,1.5*Nsave./Nsave(1),'--k','linewidth',pLwidth,'markersize',ms,'MarkerFaceColor','k'),grid on
            ax = gca;
            ax.XAxis.FontSize = ticksize;
            ax.YAxis.FontSize = ticksize;
            set(gca,'LineWidth',Lwidth)
            xlabel('\mbox{\boldmath $N$}','interpreter','latex','fontsize',lblsize)
            ylabel('\textbf{Normalized wall clock time}','interpreter','latex','fontsize',lblsize)
            % xlim([1e4 1e6])
            % xticks([1e4 1e5 1e6])
            % xticklabels({'1e4','1e5','1e6'})
            xlim([1e3 1e6])
            xticks([1e3 1e4 1e5 1e6])
            xticklabels({'1e3','1e4','1e5','1e6'})


            Nlegend =   [Nlegend;{['$N = ' num2str(Nsave(jjj)) '$']}];

        end
        subplot(1,4,1),legend(Nlegend','Location','northeast','interpreter','latex','fontsize',lgndsize),grid on
        subplot(1,4,4),legend('ML','$\mathcal{O}(N)$','Location','southeast','interpreter','latex','fontsize',lgndsize),grid on
        subplot(1,4,3),legend('ML',['$\mathcal{O}\left( \rho_{mean}^' num2str(rbf_set.poly) '\right)$'],'Location','northeast','interpreter','latex','fontsize',lgndsize),grid on

        % subplot(4,4,subpos(iii,1)),legend(Nlegend','Location','northeast','interpreter','latex','fontsize',lgndsize),grid on
        % subplot(4,4,subpos(iii,4)),legend('ML','$\mathcal{O}(N)$','Location','southeast','interpreter','latex','fontsize',lgndsize),grid on
        % subplot(4,4,subpos(iii,3)),legend('ML',['$\mathcal{O}\left( \rho_{mean}^' num2str(rbf_set.poly) '\right)$'],'Location','northeast','interpreter','latex','fontsize',lgndsize),grid on

        set(gcf,'Position',[21 230.6000 1.8688e+03 457.6000])

        text(-4.10,0.5,['$m_L$ = ' num2str(rbf_set.poly)],'interpreter','latex','units','normalized','fontsize',24)

        % set(gcf,'Position',[401.8000 199.4000 1.1064e+03 904.0000])
        drawnow
        if SAVE == 1
            filefig = ['ml2_' num2str(rbf_set.poly) '_' EQ '.fig']
            fileeps = ['ml2_' num2str(rbf_set.poly) '_' EQ '.eps']
            savefig(gcf,filefig)
            saveas(gca, fileeps,'epsc');
        end


        %polynomial iteration end


        if 1==1

            %% Error plots for all node resolutions
            figure(122)
            % tiledlayout(4,size(Xsave,2));
            % t = tiledlayout(1,size(Xsave,2),'TileSpacing','Compact');
            if iii==1
                % t = tiledlayout(4,size(Xsave,2),'TileSpacing','Compact','Padding','Compact');
                t = tiledlayout(4,size(Xsave,2),'TileSpacing','Compact');
            end
            ms = 5;
            for jjj=1:size(Xsave,2)
                % c = abs(u-uorig)./max(abs(uorig));
                c = Esave{jjj}./max(Esave{jjj});
                % subplot(4,size(Xsave,2),jjj)
                nexttile
                scatter(Xsave{jjj}(:,1),Xsave{jjj}(:,2),ms,c,'filled'),axis equal tight off,
                xlim([-0.5 0.5]),ylim([-0.5 0.5])
                % xlim([0 1]),ylim([0 1])
                if iii==1
                    title(['$N$ = ' num2str(Nsave(jjj))],'interpreter','latex','fontsize',24)
                end
                if iii==4
                    cbh = colorbar;
                    % cbh.Location = 'eastoutside';
                    cbh.Location = 'southoutside';
                    % cbh.Layout.Tile = 'northoutside';
                    cbh.Ticks = linspace(0, 1, 3);
                    cbh.Limits = [0 1];
                    cbh.FontSize = 18;
                end
                if jjj==1,ylabel('y','interpreter','latex'),end
                xlabel('x','interpreter','latex')
            end
            % cbh = colorbar ; %Create Colorbar
            set(gcf,'Position',[273 161 1624 1035])
            % set(gcf,'Position',[1 84.2000 1.4832e+03 1.0344e+03])
            % fontsize(gcf, 18, 'points')
            text(-6.74,0.5,['$m_L$ = ' num2str(rbf_set.poly)],'fontsize',24,'interpreter','latex','units','normalized' )
            % title(t,'Normalized relative errors','interpreter','latex','fontsize',30)
            drawnow
            if SAVE == 1 && iii == 4
                %     filefig = ['ml3_all' num2str(rbf_set.poly) '.fig']
                %     fileeps = ['ml3_all' num2str(rbf_set.poly) '.eps']
                savefig(gcf,['ml3_' EQ '.fig'])
                saveas(gca, ['ml3_' EQ '.eps'],'epsc');
            end

        end

    end

end

if 1==0
    %% Node coarsening animation
    figure
    set(gcf,'position',[1.8000 41.8000 766.4000 740.8000])
    PAUSE = 0;
    ms = 4;
    for rerun=1:10
        for i=1:size(L,2)
            if PAUSE == 1,pause,else,pause(.75),end
            %         figure(4)
            plot(Xlvls{i}(:,1),Xlvls{i}(:,2),'.k','MarkerSize',1.0*ms),hold on
            plot(Xlvls{i}(ismember(Xlvls{i},Xb,'rows'),1),Xlvls{i}(ismember(Xlvls{i},Xb,'rows'),2),'r.','MarkerSize',1.0*ms)
            hold off
            title(['Level: ' num2str(i)])
            axis equal tight off
            xlabel('x')
            ylabel('y')
            zlabel('level')
            drawnow
        end
    end
end









