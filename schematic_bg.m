addpath ../../
clf
plot(cumsum(randn(1,500)),'color',.8*[1 1 1],'linewidth',2)
axis off
set(gcf,'color','w')
axis tight

addpath ../../export_fig/

export_fig('-pdf','Figs/schematic_bg')