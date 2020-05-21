plot(randn(100,1)); 
set(gca,'Units','normalized','Position',[0.13 0.11 0.775 0.815]);
set(gcf,'Units','pixels','Position',[4 4 1200 900]);  %# Modify figure size
hgexport(gcf,'myfig.png',...
    hgexport('factorystyle'),'Format','png');
