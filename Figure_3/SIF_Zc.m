% This script is used to plot the stress intensity factor 
% based on finite element solution in FEniCS
% Cartesian Coordinate (x,y)
% Date: 24/02/2022
% Author: Hanwen Zhang

%% Dimensional Parameters
H_real=100; % m
rho_real=917; % kg/m^3
g_real=9.8; % m/s^2
KIc_real=100; % kPa*m^(1/2)

%% Import the value in SIF_deltatau02.txt
% when W'=1, Delta tau'=0.2, f=0.7
K_file_1 = fopen('SIF_deltatau02.txt','r');
formatSpec = '%f %f %f %f';
sizeA = [4 Inf];
A = fscanf(K_file_1,formatSpec,sizeA);
fclose(K_file_1);

vlambda_1=A(1,:);
K_Ic_1=A(2,:);
K_scale_1=K_Ic_1(1);
K1_DCM_1=A(3,:);
K2_DCM_1=A(4,:);
theta_1 = zeros(size(K1_DCM_1));
K1_Op_1 = zeros(size(K1_DCM_1));
[~,Nlambda_1]=size(K1_DCM_1);

%% Import the value in SIF_deltatau03.txt
% when W'=1, Delta tau'=0.3, f=0.7
K_file_2 = fopen('SIF_deltatau03.txt','r');
B = fscanf(K_file_2,formatSpec,sizeA);
fclose(K_file_2);
% convert it to array
vlambda_2=B(1,:);
K_Ic_2=B(2,:);
K_scale_2=K_Ic_2(1);
K1_DCM_2=B(3,:);
K2_DCM_2=B(4,:);
theta_2 = zeros(size(K1_DCM_2));
K1_Op_2 = zeros(size(K1_DCM_2));
[~,Nlambda_2]=size(K1_DCM_2);
%% theta_op and K1_op
for i=1:Nlambda_1
    theta_1(i)=-2*atan((-K1_DCM_1(i)+sqrt(K1_DCM_1(i).^2+8.*K2_DCM_1(i).^2))./(4*K2_DCM_1(i)));
    K1_Op_1(i) = cos(theta_1(i)/2)*(K1_DCM_1(i)*cos(theta_1(i)/2)^2-1.5*K2_DCM_1(i)*sin(theta_1(i)));
end
for i=1:Nlambda_2
    theta_2(i)=-2*atan((-K1_DCM_2(i)+sqrt(K1_DCM_2(i).^2+8.*K2_DCM_2(i).^2))./(4*K2_DCM_2(i)));
    K1_Op_2(i) = cos(theta_2(i)/2)*(K1_DCM_2(i)*cos(theta_2(i)/2)^2-1.5*K2_DCM_2(i)*sin(theta_2(i)));
end
%% Find K1 > 0
for i=1:Nlambda_1
    if(K1_DCM_1(i)<0)
        Nthetamax_1=i;
        break
    end
end
for i=1:Nlambda_2
    if(K1_DCM_2(i)<0)
        Nthetamax_2=i;
        break
    end
end

%% Plot K1 and K2 (normalised by Kc) as a function of crack length
figure();
FS=22;
p1=plot(K1_DCM_1(1:Nthetamax_1)/K_scale_1,vlambda_1(1:Nthetamax_1));
set(p1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineWidth',1.5,...
    'Color',[1 0 0]);
hold on
p2=plot(-K2_DCM_1(1:Nthetamax_1)/K_scale_1,vlambda_1(1:Nthetamax_1));
set(p2,'LineStyle','--',...
    'LineWidth',2,...
    'Color',[1 0 0]);
p3=plot(K1_DCM_2(1:Nthetamax_2)/K_scale_1,vlambda_2(1:Nthetamax_2));
set(p3,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineWidth',1.5,...
    'Color',[0 0 1]);
p4=plot(-K2_DCM_2(1:Nthetamax_2)/K_scale_1,vlambda_2(1:Nthetamax_2));
set(p4,'LineStyle','--',...
    'LineWidth',2,...
    'Color',[0 0 1]);
p5=plot(ones(size(0:0.02:1.0)),0:0.02:1.0);
set(p5,'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0 0 0]);

% legend
h=legend([p1,p2,p3,p4,p5],'$K_{I}^{\prime}|_{\Delta\tau^{\prime}=0.2}$',...
                        '$-K_{II}^{\prime}|_{\Delta\tau^{\prime}=0.2}$',...
                        '$K_{I}^{\prime}|_{\Delta\tau^{\prime}=0.3}$',...
                        '$-K_{II}^{\prime}|_{\Delta\tau^{\prime}=0.3}$',...
                        '$K_{I,c}^{\prime}$');
h.Interpreter='latex';
h.FontSize=FS+2;

% axis properties
ax1=gca;
ax1.FontSize=FS;
ax1.LineWidth=2;
set(gcf,'unit','normalized','position',[0.0,0.0,0.25,0.75]);
xlabel('$K/K_{Ic}$','Interpreter','latex','FontSize',FS+2);
ylabel('$Z_{C}^{\prime}$','Interpreter','latex','FontSize',FS+2,'Rotation',0,'position',[-2,0.52]);
set(gca,'XScale','line','YScale','line');
ax1.TickLabelInterpreter='latex';

xlim([0,0.2/K_scale_1]);
ylim([0,1]);
grid on
print('K_lambda_normalised','-dpng');