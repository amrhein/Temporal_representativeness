L = 50;
t0 = 21000;
tyMIN = 10;
tyMAX = 4000;
txMIN = 10;
txMAX = 4000;

tyR = linspace(tyMIN,tyMAX,L);
txR = linspace(txMIN,txMAX,L);


ty = repmat(tyR,L,1);
tx = repmat(txR,L,1)';
t0v = repmat(t0,L,L);

beta = .5;
out5 = integrate_powerlaw(ty(:),tx(:),t0v(:),beta);

beta = 1.5;
out15 = integrate_powerlaw(ty(:),tx(:),t0v(:),beta);

% Option to not plot cases where tbv>tsv
%notuse = ty(:)>tx(:);
%ty(notuse) = []; tx(notuse) = []; t0v(notuse) = [];


%%
addpath ../../export_fig
outvec = {'5' '15'};
bvec = {'0.5','1.5'};
clabv = {[0.05,.2,0.4,.6,0.8,2,5,10]
         [0.01,0.05,.1,.2,.3,0.4,0.8,2,5,10]
         %[0.01,0.1:0.2:2]
         };


for ii = 1:length(outvec)

    figure()    
    set(gcf,'color','w','position',[440   518   403   280])

 %   plout = zeros(L,L)-eps;
 %   plout(~notuse) = eval(['out' outvec{ii}]);
    plout(:) = eval(['out' outvec{ii}]);
    hold all
    pcolor(tyR,txR,log10(plout+eps)),shading interp
    caxis([-1.5,2])
    [C,h] = contour(tyR,txR,plout,clabv{ii},'k');
    colormap(flipud(hot))
    clabel(C,h,'LabelSpacing',500)
    ylabel('\tau_x (years)','fontsize',12)
    xlabel('\tau_y (years)','fontsize',12)
    grid on
    axis tight
    axis square
    plot([1000 1000],[0,4000],'k-','LineWidth',1.5)
    %text(400,200,['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    set(gca,'XTick',get(gca,'YTick'))
    set(gca,'XTickLabelRotation',45,'fontsize',12)
    set(gca,'YTick',get(gca,'XTick'))
    export_fig('-pdf',['Figs/errors_wrt_beta_' char(outvec(ii)) '_tau0_' num2str(t0)])    
    
end


