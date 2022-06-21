% This code is used to calculate the stress intensity factor 
% based on finite element solution in uncracked problem
% Cartesian Coordinate (x,y)
% frequency m
% half-width of sticky patch w
% Date: 26 11 2020

%% Dimensional Parameters
H_real=100; % m
rho_real=917; % kg/m^3
g_real=9.8; % m/s^2
KIc_real=100; % kPa*m^(1/2)

%% Import the value in SIF_deltatau02.txt
K_file = fopen('SIF_WF.txt','r');
formatSpec = '%f %f %f %f';
sizeA = [4 Inf];
A = fscanf(K_file,formatSpec,sizeA);
fclose(K_file);

vlambda= A(1,:);
K_Ic=A(2,:);
K_scale=K_Ic(1);
K1_WF = A(3,:);
K2_WF = A(4,:);
theta_WF = zeros(size(K1_WF));
K1_Op_WF = zeros(size(K1_WF));
[~,Nlambda]=size(K1_WF);

K_file_1 = fopen('SIF_DCM.txt','r');
B = fscanf(K_file_1,formatSpec,sizeA);
fclose(K_file_1);

vlambda_1 = B(1,:);
K_Ic_1=B(2,:);
K_scale_1 = K_Ic_1(1);
K1_DCM = B(3,:);
K2_DCM = B(4,:);
theta_DCM = zeros(size(K1_DCM));
K1_Op_DCM = zeros(size(K1_DCM));
[~,Nlambda_1]=size(K1_DCM);

% %% Weight Function Method
% for i=1:Nlambda
%     %% Post Processing
%     % In LEFM theory, K1 must be positive, while K2 and K3 can be negative
%     % K1_Greens(i)=max(K1_Greens(i),0);
%     %% angle of paopagation
%     theta_WF(i)=-2*atan((-K1_WF(i)+sqrt(K1_WF(i).^2+8.*K2_WF(i).^2))./(4*K2_WF(i)));
%     K1_Op_WF(i) = cos(theta_WF(i)/2)*(K1_WF(i)*cos(theta_WF(i)/2)^2-1.5*K2_WF(i)*sin(theta_WF(i)));
% end
% 
% for i=1:Nlambda_1
%     %% Post Processing
%     % In LEFM theory, K1 must be positive, while K2 and K3 can be negative
%     % K1_Greens(i)=max(K1_Greens(i),0);
%     %% angle of paopagation
%     theta_DCM(i)=-2*atan((-K1_DCM(i)+sqrt(K1_DCM(i).^2+8.*K2_DCM(i).^2))./(4*K2_DCM(i)));
%     K1_Op_DCM(i) = cos(theta_DCM(i)/2)*(K1_DCM(i)*cos(theta_DCM(i)/2)^2-1.5*K2_DCM(i)*sin(theta_DCM(i)));
% end

%% Plot theta as a function of crack length
figure();
subplot(1,2,1);
FS=18;
Markersize=8.0;
p1=plot(K1_WF,vlambda);
set(p1,'LineStyle','-',...
    'LineWidth',2.0,...
    'Color',[1 0 0]);
hold on

% p2=plot(-K2_WF(1:Ntheta),vlambda(1:Ntheta));
% set(p2,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],...
%     'MarkerSize',Markersize,...
%     'Marker','o',...
%     'LineStyle','--',...
%     'LineWidth',2.0,...
%     'Color',[1 0 0]);

p3=plot(K1_DCM,vlambda_1,'o');
set(p3,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',Markersize,...
    'Marker','o');

% p4=plot(-K2_DCM(1:Ntheta1),vlambda_1(1:Ntheta1));
% set(p4,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
%     'MarkerSize',Markersize,...
%     'Marker','o',...
%     'LineStyle','--',...
%     'LineWidth',2.0,...
%     'Color',[0 0 1]);

p5=plot(K_scale.*ones(size(0:0.02:1.0)),0:0.02:1.0);
set(p5,'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0 0 0]);



% t=title(['$\mathbf{K}$ with ','$w=$',num2str(w),', $\Delta\mu$=',...
%    num2str(deltamu),', $fl$=',num2str(fl)]);
% t.Interpreter='latex';

%h=legend([p1,p2,p3,p4],'$K_{1}$','$K_{2}$','$K_{1,optimal}$','$K_{1c}$');
h=legend([p1,p3,p5],'$K_{I}|_{WF}$',...
                        '$K_{I}|_{DCM}$',...
                        '$K_{I,c}$');
h.Interpreter='latex';
ax1=gca;
ax1.FontSize=FS;
ax1.LineWidth=2;
xlabel('$K_{I}$','Interpreter','latex','FontSize',FS+2);
ylabel('$Z_{c}$','Interpreter','latex','FontSize',FS+2,'Rotation',0);
set(gca,'XScale','line','YScale','line');
ax1.TickLabelInterpreter='latex';

% xlim([-1.5,0]);
ylim([0,1]);
grid on
print('K_lambda','-dpng');

%% Plot theta as a function of crack length
subplot(1,2,2);
FS=18;
Markersize=8.0;
p1=plot(K2_WF,vlambda);
set(p1,'LineStyle','-',...
    'LineWidth',2.0,...
    'Color',[1 0 0]);
hold on

p2=plot(K2_DCM,vlambda_1,'o');
set(p2,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',Markersize,...
    'Marker','o');


% t=title(['$\mathbf{K}$ with ','$w=$',num2str(w),', $\Delta\mu$=',...
%    num2str(deltamu),', $fl$=',num2str(fl)]);
% t.Interpreter='latex';
h=legend([p1,p2,p3],'$K_{II}|_{WF}$',...
                        '$K_{II}|_{DCM}$');
h.Interpreter='latex';
ax1=gca;
ax1.FontSize=FS;
ax1.LineWidth=2;
set(gcf,'unit','normalized','position',[0.2,0.2,0.3,0.6]);
xlabel('$K_{II}$','Interpreter','latex','FontSize',FS+2);
% ylabel('$Z_{c}$','Interpreter','latex','FontSize',FS+2,'Rotation',0);
set(gca,'XScale','line','YScale','line');
ax1.TickLabelInterpreter='latex';

% xlim([-1.5,0]);
ylim([0,1]);
grid on
print('K_lambda_benchmark','-dpng');