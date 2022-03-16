%Computes stress on grid of points: X and Y.

function [X,Y,Sxx,Syy,Sxy,Ux,Uy]=MainFrame2DFracture_StressesOnlyHNVD(H,W,L,nu,mu,noels,rhow,rhoi,g,f,X,Y,halfspaceflag,deltatau,Walls)
% % Example inputs
% H=100;        %Height of glacier
% W=H;          %Width of sticky patch
% L=10*H;       %Length of domain
% nu=0.33;      %Poisson's ratio
% mu=E/(2*(1+nu)); %Shear modulus (E is Young's)
% 
% noels=400;    %Sampling for model (discritisation)
% rhow=1e3;     %density of water
% rhoi=0.917e3; %density of ice
% g=9.8;        %accel due to grav
% f=0.7;        %flotation factor
% halfspaceflag=1; %If we want to compute in half-space
% deltatau=0.3; %sticky 'friction'
%
%
% %Creating grid with user defined sampling
% smpl=25;%100
% X = linspace(0,W*2,smpl);  
% Y = linspace(-H,0,smpl); 
% [X,Y] = meshgrid(X,Y);  %Location of observation points to compute stresses at.

% if L>W*2
%     x=[linspace(-L,-W*2,noels)];    
%     x=[x(1:end-1),linspace(-W*2,W,noels) ];
%     x2=[linspace(W,W*2,noels)];
%     x2=[x2(1:end-1),linspace(W*2,L,noels) ];
% else
%     x=[linspace(-L,0,noels)];    
%     x2=[linspace(0,L,noels)];
% end
% y=zeros(1,numel(x)); 
% y2=zeros(1,numel(x2)); 

x=linspace(-L,L,noels);
y=zeros(size(x));

sz=numel(x);
Pointsxy=[x;y]';
mystruct.line1=(1:sz);


% vertfrac=0;
% lenZH=0.3;
% if vertfrac==1
%     step=diff(x2);
%     num=numel(0:step:lenZH*H);
%     y3=linspace(0,lenZH*H,num);
%     x3=zeros(size(y3))+W;
%     PointsxyF=[fliplr(x2),x3(2:end);fliplr(y2),y3(2:end)]';
% else
%     PointsxyF=[fliplr(x2);fliplr(y2)]';
% end
% [Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
% BoundaryFlag=BoundaryFlag*0;

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% % Drawing fracture
% %Normal els as red.
% PlotFracture( P1,P2,'r' )
% %Fixed elements shown as blue
% id=find(Fdisp); 
% PlotFracture(P1(id,:),P2(id,:),'b' )
% %Titles
% title('fractures'), xlabel('x'), ylabel('y');axis equal
% %Drawing Normals
% quiver(MidPoint(:,1),MidPoint(:,2),LineNormalVector(:,1),LineNormalVector(:,2))



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option includes frictional contact
    %properties on the fault surface, elements cannot interpenetrate and
    %slip is reduced by the frictional parameters. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%

%Put to 1 to define the stresses defined in 'stress input' as strain values    
strain=0;    
    
	%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%

Pw=f*(rhoi*g*H); %Water pressure at base of glacier
Hw=Pw/(rhow*g); %Rearragne eq.7;

z=MidPoint(:,2);

Hw=Pw/(rhow*g); %Rearrange eq.7;
depthofwatertop=H-Hw;
%Measured from surface
zIce=z;%linspace(-H,0,100);
%Measure from water top
zWater=z;%linspace(-H,0,100);
zWater=zWater+depthofwatertop;
zWater(zWater>0)=0;
%Neg as 'z' is measured from surface
waterpressure=-(zWater*rhow*g);
iceweight=(zIce*rhoi*g);
pressure=waterpressure+iceweight;

Sxx = zeros(size(MidPoint(:,1))); 				
Syy = Sxx;%ones(size(MidPoint(:,1))).*(rhoi*g*MidPoint(:,2));
StickyPatch=abs(MidPoint(:,1))<W;
%Forces
deltashr=deltatau*(rhoi*g*H);
fexcess=-deltashr;
fbase=+W/L*deltashr;

Sxy=Sxx;
Sxy(~StickyPatch) = fbase;      
Sxy(StickyPatch) = fexcess+fbase;  

% FrakFlag=z>-H & abs(MidPoint(:,1))<L;
% Sxx(FrakFlag)=pressure(FrakFlag); %Pw
% Syy(FrakFlag)=pressure(FrakFlag);
% Sxy(FrakFlag)=0;

if Walls==1
WallFlag=abs(MidPoint(:,1))==L;
Sxx(WallFlag)=0; %Pw
Syy(WallFlag)=0;
Sxy(WallFlag)=0;
end

% Frictional parameters, define as single value for all elements or vary
% based on element number
Mu  = 0;       Mu=repmat(Mu,numel(MidPoint(:,1)),1);  %Coefficient of friction
Sf  = 0;       Sf=repmat(Sf,numel(MidPoint(:,1)),1);  %Frictional strength
%  
% figure;
% scatter(MidPoint(:,1),MidPoint(:,2),15,Sxy,'filled');

Option='C'; 

bottomboundflag=MidPoint(:,2)==-H;
Fdisp=bottomboundflag;
Option='H'; %No displacement along bottom boundary (Fdisp flag) 


%Flag to say we min strain energy release rate
StrainEnergyDir=0;
ThetaIn=NaN;  %The list of directions we are going to check (degrees from the current tip angle)
noNewels=NaN;
Indx=NaN;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(Option,'A')
    [Ds,Dn,Sxx,Syy,Sxy,~]=SlipCalculator2dXtraCond(MidPoint,HalfLength,Sxx,Syy,Sxy...
        ,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option,...
        P1,P2,StrainEnergyDir,ThetaIn,Indx,noNewels);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% % Plot displacement discontinuity magnitude for all elements
% GraphSlipVsElNum( Dn,Ds )
% 
% % Creating directions and magnitudes of slip for plotting
% PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )

 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% X=MidPoint(:,1);
% Y=MidPoint(:,2);
[dimx,dimy] = size(X);  
X=X(:);
Y=Y(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dist=0.05;
% [X,Y]=NanTolDistLine2Pnt( X,Y,P1,P2,LineNormalVector,Dist );
% %Drawing these and the location compared to the Obs Points
% figure;
% PlotFracture( P1,P2,'r' )
% title('Fractures and Observation points');
% xlabel('x'), ylabel('y');axis equal
% hold on
% scatter(X(:),Y(:),'.')
% hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
= CalculateStressesOnSurroundingPoints2d(X,Y,MidPoint,HalfLength...
,Sxx,Syy,Sxy,nu,E,halfspace,Ds,Dn,LineNormalVector);
 
%Extracting stress
[ Sxx,Syy,Sxy ] = ExtractCols( StressTChg );


%Reshaping stresses to grid dimensions
[X,Y,Sxx,Syy,Sxy]=ReshapeData2d...
    ( dimx,dimy,X,Y,Sxx,Syy,Sxy );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d...
(X,Y,MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);

end
