%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%Conventions:
%
% This follows the same unit convention as Poly3d so ''assumes
% dimensionally consistent units for physical quantities with
% dimensions of length or stress. Thus if you use kilometers for any
% quantity with dimensions of length (e.g. a coordinate position), you
% must use kilometers for all quantities with dimensions of length
% (e.g. displacement discontinuity).'' Andrew Lyell Thomas, Masters
% thesis, Stanford University.
%
% Engineering stress convention is used in this script for both tensors
% and principal stresses. S1 is the LEAST compressive stress. 
% See figure 6.13 & 6.14 in Pollard and Fletcher 2005. 

%Clearing old figures and arrays. 
clear;close all

%Loading a colourmap to be used for figures.
cmap = colormap_cpt('Ccool-warm');
cmap2 = colormap_cpt('Ccool-warm2');

     
%% Inputs and lengths (of model):
H=1000;              % Domain height
W=100;              % Sticky width
L=10*W;             % Domain length
nu=0.33;            % Poisson's ratio
E=10e9;             % Young's mod
mu=E/(2*(1+nu));    % Shear mo
Kc=0.1e6;           % Fracture toughness
rhow=1e3;           % Density of water
rhoi=0.917e3;       % Density of ice
g=9.8;              % Accel due to grav
f=0.7;              % Flotation factor
deltatau=0.2;       % DeltaTau - 'sticky friction'
halfspaceflag=1;    % Compute in a half-space?
Walls=1;            % Flag to turn on the lateral 'traction-free' walls   

%% Numerical bits (sampling, method and start geom):
noels=1250;         % Relates to the no of BEM elements used in calc
lenZH=0.03;         % Length of the initial fracture (as fraction of H): 
                    % The fracture starts vertical, 
                    % if you want this shorter decrease this and increase
                    % the sampling 'no els' -> note longer simulation
                    % times. I have added an error if this is less than 5
                    % els long.
StrainEnergyDir=1;  % Flag to say we use min strain energy release rate or 
                    % K1 K2 ratios to compute next propagation dir.
grdsmpl=45;         % Sampling of obs points for simple calc.

%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(-L,L,noels);
y=zeros(size(x));
%Make sure not cutting through el midpoint with fracture tip.
[val,ToMv]=min(abs(x-W));
step=diff(x);
x=x+val;%-step(1)/2;
x(end)=L;
x(1)=-L;
grad=gradient(x);
if grad(1)<grad(2)/3
    error('try different sampling')
end
if grad(end)<grad(2)/3
    error('try different sampling')
end
sz=numel(x);
Pointsxy=[x;y]';
mystruct.line1=(1:sz);


vertfrac=1;
if vertfrac==1
    step=diff(x);
    num=numel(0:step(1)/3:lenZH*W);
    y3=linspace(0,lenZH*W,num);
    x3=zeros(size(y3))+W;
    if numel(y3)<5
        disp(y3)
        %error('You need a higher sampling')
    end
    PointsxyF=[x3(1:end);y3(1:end)]';
end
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
BoundaryFlag=BoundaryFlag*0;


if Walls==1
    %Add walls:
    step=diff(x);
    num=numel(0:step:H);
    y4=linspace(0,H,num);
    x4=zeros(size(y4))+L;
    PointsxyF=[(x4);(y4)]';
    [Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
    
    x4=zeros(size(y4))-L;
    PointsxyF=[fliplr(x4);fliplr(y4)]';
    [Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
    BoundaryFlag=BoundaryFlag*0;
end


%Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
    Fdisp=[];
end
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = halfspaceflag; 

% Defining the height of the freesurface relative to the value of 0 of
% your imported surfaces. 
freesurface_height = H;
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;


% Defining elastic constants
% mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
% nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Reshaping the list of XY into usable variables
[MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp  ]...
 = CreateElements2d( Pointsxy,mystruct,Fdisp );

% Drawing fracture
%Normal els as red.
PlotFracture( P1,P2,'r' )
%Fixed elements shown as blue
id=find(Fdisp); 
PlotFracture(P1(id,:),P2(id,:),'b' )
%Titles
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Drawing Normals
quiver(MidPoint(:,1),MidPoint(:,2),LineNormalVector(:,1),LineNormalVector(:,2))
scatter(x,ones(size(x))*-H,'k')
no=100;
KtotalLp=zeros(no,1);
KILp=zeros(no,1);
K2Lp=zeros(no,1);
GLp=zeros(no,1);
EndPointLp=zeros(no,2);

%% Now compute stresses without fracture and path based on most compressive stress
%direction:
X = linspace(0,2*W,grdsmpl);  
Y = linspace(-H+100*lenZH,-H+W,grdsmpl); 
[X,Y] = meshgrid(X,Y);  %Location of observation points to compute stresses at.
[dimx,dimy] = size(X);  
[X,Y,Sxx,Syy,Sxy,Ux,Uy]=MainFrame2DFracture_StressesOnlyHNVD(H,W,L,nu,mu,noels,rhow,rhoi,g,f,X,Y,halfspaceflag,deltatau,Walls);


Pw=f*(rhoi*g*H); %Water pressure at base of glacier
%Water pressure/lithostatic on fracture
Hw=Pw/(rhow*g); %Rearrange eq.7;
depthofwatertop=H-Hw;
%Measured from surface
zIce=Y;%linspace(-H,0,100);
%Measure from water top
zWater=Y;%linspace(-H,0,100);
zWater=zWater+depthofwatertop;
zWater(zWater>0)=0;
%Neg as 'z' is measured from surface
waterpressure=-(zWater*rhow*g);
iceweight=(zIce*rhoi*g);
pressure=waterpressure+iceweight;

Sxx=Sxx;
Syy=Syy;

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

%Drawing S1 S2 directions.
%DrawS1S2Directions(X(:),Y(:),S1dir,S2dir );hold on
S2x=S2dir(:,1);
S2y=S2dir(:,2);
%Flip
S2x(S2y<0)=-S2x(S2y<0);
S2y(S2y<0)=-S2y(S2y<0);

[X,Y,S2x,S2y]=ReshapeData2d...
    ( dimx,dimy,X,Y,S2x,S2y );

save('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=1000_W=100_L=1000_f=0_7.mat',...
    'X','Y','S1','S2','S2x','S2y','W','H','lenZH','-append');
save('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_2_H=1000_W=100_L=1000_f=0_7.mat',...
    'X','Y','S1','S2','S2x','S2y','W','H','lenZH','-append');
save('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_1_H=1000_W=100_L=1000_f=0_9.mat',...
    'X','Y','S1','S2','S2x','S2y','W','H','lenZH','-append');