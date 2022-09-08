   

    function plot_Peanut_withboundary(filename,savename,minimum,maximum)
        plotarray = importdata(filename);
        x = linspace(0,1,301);
        y = x;
        n=100;
        ii = 0:1:n;
        ti = 2*pi*ii/n;
        xfun = @(x) 0.06*((cos(x)+2).*(cos(x+0.6)+2).*(2+0.1*(cos(3*x))))-0.1;
        yfun = @(y) 0.06*((sin(y)+2).*(sin(y-0.5)+2).*(2+0.4*(cos(2*y))).*(1+0.1*(sin(4*y))))-0.06;
        vix = xfun(ti);
        viy = yfun(ti);
        figure(1)
        hold on
        Z = real(plotarray);
        imagesc(x,y,Z);
        ax = gca;
        set(gca,'Yscale','linear')
        set(gca,'Xscale','linear')
        set(gca,'visible','off')
        set(ax, 'XLimMode', 'auto', 'YLimMode', 'auto')
        plot(vix,viy,'k','LineWidth',2.5);
        axis equal
        axis([0 1 0 1])
        colorbar;
        set(gca,'YDir','normal')
        caxis( [minimum-0.00001 maximum] )
        cmap = [1 1 1; parula(10000)];
        colormap(cmap);
        saveas(gcf,savename,'epsc')
        hold off
    end
    %end
