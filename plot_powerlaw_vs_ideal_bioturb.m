L = 50;
t0 = 20000;
td = 400;
tbMIN = 10;
tbMAX = 2000;
tsMIN = 10;
tsMAX = 2000;

tbR = linspace(tbMIN,tbMAX,L);
tsR = linspace(tsMIN,tsMAX,L);

tbv = repmat(tbR,L,1);
tsv = repmat(tsR,L,1)';
t0v = repmat(t0,L,L);
tdv = repmat(td,L,L);

bet = .5;
out5 = integrate_powerlaw_vs_ideal_bioturb(tbv,tsv,tdv,t0v,bet);

bet = 1.5;
out15 = integrate_powerlaw_vs_ideal_bioturb(tbv,tsv,tdv,t0v,bet);

%%
addpath ../../export_fig
close all
outvec = {'5' '15'};
bvec = {'0.5','1.5'};

for ii = 1:length(outvec)

    figure(ii)
    set(gcf,'color','w','position',[440   518   403   280])
    
    ploutv = eval(['out' outvec{ii}]);
    plout = reshape(ploutv,length(tbR),length(tsR));
    
    hold on

    %mmi = find(plout(:,1)==min(plout(:,1)));
    %mma = find(plout(:,end)==min(plout(:,end)));
    %[~,mm] = min(plout)
    % Plots the minimum values for each ts as a line
    %plot(tbR,tsR(mm),'-k','linewidth',2)
    %plot([min(tbR),max(tbR)],[tbR(mmi),tsR(mma)],'--k','linewidth',2)
    %plot([min(tbR),max(tbR)],[min(tsR),max(tsR)],'-k','linewidth',2)
    
    % Plots the minimum values for each ts as a matrix
    %imagescnan(tbR,tsR,double(bsxfun(@eq,plout,min(plout))))
    
    % All of plout as a matrix
    %imagescnan(tbR,tsR,(plout));
    %h = colorbar;
    %caxis([0,.5])
    %colormap(flipud(gray))

    [C,h] = contour(tbR,tsR,plout,[0.02,0.05,0.1:0.1:.8],'k');
    %contour(tbR,(tsR(2:end)+tsR(1:end-1))/2,diff(plout),[0,0],'k','linewidth',2);
    contour(tbR,tsR(2:end),diff(plout),[0,0],'k','linewidth',2);
    clabel(C,h)
    grid on
    ylabel('\tau_s (years)','fontsize',12)
    xlabel('\tau_b (years)','fontsize',12)
    axis tight
    axis square
    %text(400,200,['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    %title(['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    set(gca,'XTick',get(gca,'YTick'))
    set(gca,'XTickLabelRotation',45,'fontsize',12)
    set(gca,'YTick',get(gca,'XTick'))
    export_fig('-pdf',['Figs/errors_wrt_beta_vs_ideal_bioturb_' char(outvec(ii)) '_tau0_' num2str(t0)])    
end

