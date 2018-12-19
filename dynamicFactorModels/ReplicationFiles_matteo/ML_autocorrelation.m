% ML_autocorrelation - Computate Auto-Correlation up to order g
%
% [acf pacf]=ML_autocorrelation(Y,g,graph);
% Inputs:
%   Y variable
%   g lag up to which compute the autocorrelation
%   plot (optional) if 1 plots autocorrelations
%

% Written by Matteo Luciani

function [acf, pacf]=ML_autocorrelation(Y,g,graph,label,nopacf)

[t, n] = size(Y);
if nargin==1; g = round(t^.5); end
if nargin<5; nopacf=0; end
acf=zeros(g,n); pacf=zeros(g,n);
for ii=1:n;
    m1=mean(Y(:,ii));
    y = Y(:,ii) - m1;
    s2 = y'*y;
    for j=1:g;
        x1=y(j+1:t);
        x2=y(1:t-j);
        acf(j,ii)=(x1'*x2)/s2;
        clear x1 x2
    end;
    clear y m1
end;

if nopacf==0;
    for jj=1:n;
        for ii=1:g;
            beta=ML_autoregressive(Y(:,jj),1,g,ii);
            pacf(ii,jj)=beta(ii+1,1);
        end
    end;
else
    pacf=zeros(1,n);
end
if g>1;
    acf=[ones(1,n); acf];   
    pacf=[ones(1,n); pacf];
end

if nargin>2;
%     stderrs=((1+2*cumsum(acf.^2))/t).^(.5); stderrs=[zeros(1,n); stderrs]; up = 1.96*stderrs; down = -1.96*stderrs;   
    up2=2/sqrt(t)*ones(g+1,1); down2=-2/sqrt(t)*ones(g+1,1); % hamilton (1994) pag 111
    tt=0:1:g;
    kk=ceil(sqrt(n)); if kk^2-n>=kk;kk2=kk-1;else kk2=kk; end   
    if graph==1;axes('Parent',figure,'FontSize',10); ML_FigureSize                       
        if n==1; 
            if g>25; 
                hold on; 
                plot(tt,up2(:,1),'r-',tt,down2(:,1),'r-');
                p1=stem(tt,[acf pacf],'filled','linewidth',1.5);   
                set(p1(1),'Color','k')
                set(p1(2),'Color',[.5 .5 .5])
                hold off;  box on
            else
                p1=bar(tt,[acf pacf]);
                set(p1(1),'faceColor','k')
                set(p1(2),'faceColor',[.5 .5 .5])
                hold on; plot(tt,up2(:,1),'r-',tt,down2(:,1),'r-'); hold off;
            end; 
            legend(p1,{'acf','pacf'});                  
            axis([0 g -1 1]);
            gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1)
        else            
            for i=1:n; 
                subplot(kk,kk2,i),
                if g>25; 
                    plot(tt,[acf(:,i) pacf(:,i)],'linewidth',1.5); 
                else
                    bar(tt,[acf(:,i) pacf(:,i)]); 
                end;  
                legend('acf','pacf','location','southeast'); axis([0 g -1 1]);             
                if nargin==4; 
                    if iscell(label);
                        title(label{i},'fontweight','bold','fontsize',12); 
                    else
                        title(label(i,:),'fontweight','bold','fontsize',12);
                    end;            
                end
                hold on; plot(tt,up2,'r-',tt,down2,'r-'); hold off; 
            end
        end
    elseif graph==2;
        for i=1:n; axes('Parent',figure,'FontSize',10); ML_FigureSize
            subplot(2,1,1),plot(Y(:,i),'linewidth',1.5); axis tight; if nargin==4; title(label(i,:),'fontweight','bold','fontsize',12); end;
            subplot(2,1,2),bar(tt,[acf(:,i) pacf(:,i)]);  legend('acf','pacf'); axis([0 g -1 1]);             
            hold on; plot(tt,up2,'r-',tt,down2,'r-'); hold off;
%             print('-dpdf','-painters','-r600',['Fattore_' num2str(i) '.pdf'])
        end;
    end
end

