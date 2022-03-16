figure(4);hold on
set(gcf,'Position',[681   557   480*2   307.5*1.5])
colorpath=['/Users/eart0487/Documents/Oxford/Oxford Research/',...
    'ice sheets/FE Model/redblue'];
addpath(colorpath);
%% Pabel(a)
% H=1000m, W=100m, for the nondimensional case, H=1,W'=0.1
ax1=subplot(1,2,1);
% BEM Simulations
% curve 1
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=1000_W=100_L=1000_f=0_7_nd.mat');
% analytical stress trajectories
lineobj =streamline(X/H,Y/H,S2x,S2y,W/H,-1+0.03*100/H,[0.1 110]);
% x=lineobj.XData;
% y=lineobj.YData;
set(lineobj,'LineWidth',2.5)
set(lineobj,'LineStyle','--')
set(lineobj,'Color','black');
% Drawing fracture
%Normal els as green.
LineType='-';
LineColor=[237/256 27/256 58/256]; %red
LineWidth=5.0;
frac11=PPlotFracture(P1/H,P2/H,LineType, LineColor,LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
PPlotFracture(P1(id,:)/H,P2(id,:)/H, LineType, 'w', LineWidth)

% curve 2
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_1_H=1000_W=100_L=1000_f=0_9_nd.mat');
% Drawing fracture
LineType='-';
LineColor=[0 201/256 87/256];
LineWidth=3.0;
frac12=PPlotFracture(P1/H,P2/H,LineType, LineColor, LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
PPlotFracture(P1(id,:)/H,P2(id,:)/H, LineType, 'w', LineWidth)

% curve 3
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_2_H=1000_W=100_L=1000_f=0_7_nd.mat');
% Drawing fracture
LineType='-';
LineColor='b'; % blue
LineWidth=2.0;
frac13=PPlotFracture(P1/H,P2/H,LineType, LineColor,LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
bottom1=PPlotFracture(P1(id,:)/H,P2(id,:)/H,'--',[0 0 0],LineWidth);

% Tick and label
FS=20;
ax1.FontSize=FS-4;
ax1.LineWidth=1.5;
ax1.TickLabelInterpreter='latex';
% axis range
set(ax1,'XScale','line','YScale','line');
axis('equal')
grdsmpl=100;
yheight=75;

% %Remove some vals of stream line?
LineLen=0.5;
%figure;hold on
S_div=S1-S2;
Sdiv_max=max(max(abs(S_div)));
quiver(X(:)/H,Y(:)/H,abs(S_div(:))./Sdiv_max.*S2x(:),abs(S_div(:))./Sdiv_max.*S2y(:),...
    LineLen,'k.')
 
txt3 = '$\mathbf{\left(a\right)}$';
t3 = text(W/H-180/H,-1+700/H,txt3,'Interpreter','latex','Fontsize',FS-2);

% WhiteFigure;
axis('equal')
xl=ylabel('$x^{\prime}$','Interpreter','latex','FontSize',FS-2,'Rotation',0);
xl.Position=[125/H -H-0.1 -1.0000];
yl=ylabel('$z^{\prime}$','Interpreter','latex','FontSize',FS-2,'Rotation',0);
yl.Position=[-175/H -1+370/H -1.0000];

xlim([W/H-0.2 W/H+0.2])
ylim([-1 -1+yheight*10/H])
yticks(-1:100/H:-1+yheight*10/H);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'})
caxis([0 25])
box on
r = rectangle('Position',[80/H -1 40/H yheight/H]');
r.LineWidth = 2.0;

set(ax1,'Position',[0.0500 0.1000 0.3000 0.8000]);

h=legend([frac11 frac13 frac12 bottom1], '$f=0.7,\ \Delta\tau^{\prime}=0.3$',...
                        '$f=0.7,\ \Delta\tau^{\prime}=0.2$',...
                        '$f=0.9,\ \Delta\tau^{\prime}=0.1$',...
                        '$\sigma_{2}\ trajectory$');
h.Interpreter='latex';
h.FontSize=FS-8;
%%
xa=0.25;
ax2=axes('position',[0.14 0.3 xa/1.3875 7.5/4.5*xa]);
box on % put box around new pair of axes
% curve 1
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=1000_W=100_L=1000_f=0_7.mat');
% analytical stress trajectories
lineobj =streamline(X/H,Y/H,S2x,S2y,W/H,-1+0.03*100/H);
x=lineobj.XData;
y=lineobj.YData;
set(lineobj,'LineWidth',2.5)
set(lineobj,'LineStyle','--')
set(lineobj,'Color','black');
% Drawing fracture
%Normal els as red.
LineType='-';
LineColor=[237/256 27/256 58/256]; %red
LineWidth=4.0;
frac21=PPlotFracture(P1/H,P2/H,LineType, LineColor,LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
PPlotFracture(P1(id,:)/H,P2(id,:)/H, LineType, 'w', LineWidth)

% curve 2
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_1_H=1000_W=100_L=1000_f=0_9.mat');
% Drawing fracture
LineType='-';
LineColor=[0 201/256 87/256];
LineWidth=2.5;
frac22=PPlotFracture(P1/H,P2/H,LineType, LineColor, LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
PPlotFracture(P1(id,:)/H,P2(id,:)/H, LineType, 'w', LineWidth)

% curve 3
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_2_H=1000_W=100_L=1000_f=0_7.mat');
% Drawing fracture
LineType='-';
LineColor='b'; % blue
LineWidth=2.0;
frac23=PPlotFracture(P1/H,P2/H,LineType, LineColor,LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
bottom2=PPlotFracture(P1(id,:)/H,P2(id,:)/H,'--',[0 0 0],LineWidth);

% Tick and label
FS=16;
ax2.FontSize=FS-4;
ax2.LineWidth=1.5;
ax2.TickLabelInterpreter='latex';
% axis range
set(ax2,'XScale','line','YScale','line');
axis('equal')
grdsmpl=100;
yheight=75;

LineLen=0.5;
%figure;hold on
S_div=S1-S2;
Sdiv_max=max(max(abs(S_div)));
quiver(X(:)/H,Y(:)/H,abs(S_div(:))./Sdiv_max.*S2x(:),abs(S_div(:))./Sdiv_max.*S2y(:),...
    LineLen,'k.')

% WhiteFigure;
axis('equal')
% yl=ylabel('$z^{\prime}$','Interpreter','latex','FontSize',FS-2,'Rotation',0);
% yl.Position=[-150/H -1+370/H -1.0000];
% xl=xlabel('$x^{\prime}$','Interpreter','latex','FontSize',FS-2);

xlim([80/H 120/H])
xticks(0.08:0.02:0.12);
xticklabels({'0.08','0.10','0.12'})
ylim([-1 -1+yheight/H])
yticks(-1:20/H:-1+yheight/H);
yticklabels({'0','0.02','0.04','0.06','0.08','0.10'})
box on
% h=legend([frac1 frac3 frac2 bottom], '$f=0.7,\ \Delta\tau^{\prime}=0.3$',...
%                         '$f=0.7,\ \Delta\tau^{\prime}=0.2$',...
%                         '$f=0.9,\ \Delta\tau^{\prime}=0.1$',...
%                         '$\sigma_{2}\ trajectory$');
% h.Interpreter='latex';
% h.FontSize=FS-8;

%% Pabel(b)
ax3=subplot(1,2,2);
yheight=0.75;

% BEM Simulations
% curve 1
% load('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=100_W=100_L=1000_f=0_7.mat');
% load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=1000_W=1000_L=10000_f=0_7.mat');
load('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=1000_W=1000_L=10000_f=0_7.mat');
% analytical stress trajectories

lineobj =streamline(X/H,Y/H,S2x,S2y,W/H,-1+0.03*1000/H);
set(lineobj,'LineWidth',2.5)
set(lineobj,'LineStyle','--')
set(lineobj,'Color','black');
% Drawing fracture
%Normal els as red.
LineType='-';
LineColor=[237/256 27/256 58/256]; %red
LineWidth=3.5;
frac31=PPlotFracture(P1/H,P2/H,LineType, LineColor,LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
PPlotFracture(P1(id,:)/H,P2(id,:)/H, LineType, 'w', LineWidth)

% curve 2
% load('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_01_Dt=0_3_H=1000_W=1000_L=10000_f=0_7.mat');
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_1_H=1000_W=1000_L=10000_f=0_9.mat');
% Drawing fracture
LineType='-';
LineColor=[0 201/256 87/256];
LineWidth=2.5;
frac32=PPlotFracture(P1/H,P2/H,LineType, LineColor, LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
PPlotFracture(P1(id,:)/H,P2(id,:)/H, LineType, 'w', LineWidth)

% curve 3
% load('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_2_H=100_W=100_L=1000_f=0_7.mat');
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_2_H=1000_W=1000_L=10000_f=0_7.mat');
% Drawing fracture
LineType='-';
LineColor='b'; % blue
LineWidth=2.0;
frac33=PPlotFracture(P1/H,P2/H,LineType, LineColor,LineWidth);
%Fixed elements shown as blue
id=find(Fdisp); 
bottom3=PPlotFracture(P1(id,:)/H,P2(id,:)/H,'--',[0 0 0],LineWidth);

% Tick and label
FS=20;
ax3.FontSize=FS-4;
ax3.LineWidth=1.5;
ax3.TickLabelInterpreter='latex';
% axis range
set(ax3,'XScale','line','YScale','line');
axis('equal')


% %Remove some vals of stream line?
LineLen=0.5;
%figure;hold on
S_div=S1-S2;
Sdiv_max=max(max(abs(S_div)));
quiver(X(:)/H,Y(:)/H,abs(S_div(:))./Sdiv_max.*S2x(:),abs(S_div(:))./Sdiv_max.*S2y(:),...
    LineLen,'k.')
 
txt1 = '$\mathbf{\left(b\right)}$';
t1 = text(W/H-0.18,-1+0.7,txt1,'Interpreter','latex','Fontsize',FS-2);

% h=legend([frac1 frac3 frac2 bottom], '$f=0.7,\ \Delta\tau^{\prime}=0.3$',...
%                         '$f=0.7,\ \Delta\tau^{\prime}=0.2$',...
%                         '$f=0.9,\ \Delta\tau^{\prime}=0.1$',...
%                         '$\sigma_{2}\ trajectory$');
% h.Interpreter='latex';
% h.FontSize=FS-8;

% WhiteFigure;
axis('equal')
xl=xlabel('$x^{\prime}$','Interpreter','latex','FontSize',FS-2,'Rotation',0);

xlim([W/H-0.2 W/H+0.2])
ylim([-1 -1+yheight])
yticks(-H/H:0.1:-1+yheight);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'})
set(ax3,'Position',[0.3500 0.100 0.3000 0.8000]);
box on

%% Pabel(c)
load('Results/NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_2_H=10_W=100_L=1000_f=0_7.mat');
ax4=subplot(1,3,3);
yheight=75;
% BEM Simulations
% curve 1
%                      这里load uncracked stress field
% analytical stress trajectories
lineobj =streamline(X/H,Y/H,S2x,S2y,W/H,-1+0.03*100/H);
set(lineobj,'LineWidth',2.5)
set(lineobj,'LineStyle','--')
set(lineobj,'Color','black');
hold on
% Vertical fracture
Z_C_1=0.559;
LineType='-';
LineColor=[237/256 27/256 58/256]; %red
LineWidth=3.5;
frac41=plot([W/H W/H],[0-1 Z_C_1-1],'LineStyle',LineType,...
    'Color',LineColor,...
    'LineWidth',LineWidth);
hold on
% Drawing fracture
Z_C_2=0.443;
LineType='-';
LineColor=[0 201/256 87/256];
LineWidth=2.5;
frac42=plot([W/H W/H],[0-1 Z_C_2-1],'LineStyle',LineType,...
    'Color',LineColor,...
    'LineWidth',LineWidth);

% Drawing fracture
Z_C_3=0.26;
LineType='-';
LineColor=[0 0 1]; % blue
LineWidth=2.0;
frac43=plot([W/H W/H],[0-1 Z_C_3-1],'LineStyle',LineType,...
    'Color',LineColor,...
    'LineWidth',LineWidth);
% Tick and label
FS=20;
ax4.FontSize=FS-4;
ax4.LineWidth=1.5;
ax4.TickLabelInterpreter='latex';
% axis range
set(ax4,'XScale','line','YScale','line');
axis('equal')

% id=find(Fdisp); 
% bottom=PPlotFracture(P1(id,:)/H,P2(id,:)/H,'--',[0 0 0],LineWidth);

% %Remove some vals of stream line?
LineLen=0.25;
%figure;hold on
S_div=S1-S2;
Sdiv_max=max(max(abs(S_div)));
quiver(X(:)/H,Y(:)/H,abs(S_div(:))./Sdiv_max.*S2x(:),abs(S_div(:))./Sdiv_max.*S2y(:),...
    LineLen,'k.')
 
txt1 = '$\mathbf{\left(c\right)}$';
t1 = text(W/H-0.18*H/H,-1+0.7,txt1,'Interpreter','latex','Fontsize',FS-2);

h=legend([frac41 frac42 frac43], '$f=0.7,\ \Delta\tau^{\prime}=0.036$',...
                        '$f=0.7,\ \Delta\tau^{\prime}=0.035$',...
                        '$f=0.7,\ \Delta\tau^{\prime}=0.033$');
h.Interpreter='latex';
h.FontSize=FS-8;

xlim([W/H-0.2 W/H+0.2])
ylim([-1 -1+yheight/H])
yticks(-H/H:10/H:-1+yheight/H);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'})
set(ax4,'Position',[0.6500 0.1000 0.3000 0.8000]);


saveas(gca,'FracturePath_BEM','jpeg')
saveas(gca,'FracturePath_BEM','fig')