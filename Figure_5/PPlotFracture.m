function pl=PPlotFracture(P1,P2,LineType, LineColor, LineWidth)
%PlotFracture Plots the loaded 2d fractures and the normals of each
%element. 

hold on    
pl=line([P1(:,1)';P2(:,1)'],[P1(:,2)';P2(:,2)'],'LineStyle',...
    LineType,'LineWidth',LineWidth, 'Color',LineColor);
pl=pl(1);

end