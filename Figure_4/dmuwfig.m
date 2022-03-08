% The code is used to visuize the results
% plot crack length d as a function of width and delta tau
% 17/12/2020
% Hanwen Zhang
vw = 0.1:0.1:1.0;
vtau = 0.1:0.05:0.5;
[Mtau,Mw]=meshgrid(vtau,vw);% mesh

% Interpolation
vw_interp = 0.1:0.02:0.9;
vtau_interp = 0.1:0.01:0.4;
[xq,yq] = meshgrid(vw_interp, vtau_interp);

% import the cracklength
M1 = readmatrix('Zc_f=0_7.txt');
% vq = griddata(Mw,Mtau,M,xq,yq);
[Ncase1,~] = size(M1);
for i=1:Ncase1
    cracklength1(round(M1(i,1))+1,round(M1(i,2))+1)= (M1(i,3));
end

% import the cracklength
M2 = readmatrix('Zc_f=0_9.txt');
% vq = griddata(Mw,Mtau,M,xq,yq);
[Ncase2,~] = size(M2);
for i=1:Ncase2
    cracklength2(round(M2(i,1))+1,round(M2(i,2))+1)= (M2(i,3));
end


%% Pcolor

figure()
subplot(1,2,1);
FS=18;
p=pcolor((Mw),(Mtau),(cracklength1));
shading faceted
hold on;
txt = '$\mathbf{\left(a\right)}$';
t = text(0.1,0.475,txt,'Interpreter','latex','Fontsize',20,'color',[1 1 1]);

pos=axis;
ylabel('$\Delta\tau^{\prime}$','Interpreter','latex','Rotation',0,'position',[-0.7*pos(1) 0.55*pos(4)]);
xlabel('$W^{\prime}$','Interpreter','latex','Rotation',0);

ax=gca;
set(gca,'XScale','linear','YScale','linear');
ax.TickLabelInterpreter='latex';
ax.FontSize=FS;
ax.LineWidth=1.5;
set(gca, 'Units', 'normalized', 'Position', [0.12 0.2 0.3 0.6])
% ylim([0,0.6]);
% xlim([0,2]);


subplot(1,2,2);
p=pcolor((Mw),(Mtau),(cracklength2));
shading faceted
hold on;
xlabel('$W^{\prime}$','Interpreter','latex','Rotation',0);

% colormap();
cbar=colorbar;
cbar.TickLabelInterpreter = 'latex';
caxis([0,1.0]);
cbar.Ticks = linspace(0,1.0,6) ; %Create 8 ticks from zero to 1
cbar.TickLabels = {'$\le 0.05$','$0.2$','$0.4$','$0.6$','$0.8$','$1.0$'} ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
ylabel(cbar,'$Z_{C,max}^{\prime}$','Interpreter','latex','Rotation',0);

ax=gca;
set(gca,'XScale','linear','YScale','linear');
ax.TickLabelInterpreter='latex';
ax.FontSize=FS;
ax.LineWidth=1.5;
% ylim([0,0.6]);
% xlim([0,2]);
set(gca, 'Units', 'normalized', 'Position', [0.5 0.2 0.3 0.6])
set(gcf,'position',[0 0 800 400]);
txt = '$\mathbf{\left(b\right)}$';
t = text(0.1,0.475,txt,'Interpreter','latex','Fontsize',20,'color',[1 1 1]);

