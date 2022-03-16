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
H=1e3;              % Domain height
W=1e3;              % Sticky width
L=10*W;             % Domain length
nu=0.33;            % Poisson's ratio
E=10e9;             % Young's mod
mu=E/(2*(1+nu));    % Shear mo
Kc=0.1e6;           % Fracture toughness
rhow=1e3;           % Density of water
rhoi=0.917e3;       % Density of ice
g=9.8;              % Accel due to grav
f=0.7;              % Flotation factor
deltatau=0.3;       % DeltaTau - 'sticky friction'
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for i=1:no

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
waterpressure=-(zWater.*rhow.*g);
iceweight=(zIce*rhoi*g);
pressure=waterpressure+iceweight;



Sxx = zeros(size(MidPoint(:,1))); 				
Syy = Sxx;%ones(size(MidPoint(:,1))).*(rhoi*g*MidPoint(:,2));
StickyPatch=abs(MidPoint(:,1))<W;

%Forces
deltashr=deltatau*(rhoi*g*H);%0.3*(rhoi*g*H);
fexcess=-deltashr;
fbase=+W/L*deltashr;

Sxy=Sxx;
Sxy(~StickyPatch) = fbase;      
Sxy(StickyPatch) = fexcess+fbase;  

FrakFlag=z>-H & abs(MidPoint(:,1))<L;
Sxx(FrakFlag)=pressure(FrakFlag); %Pw
Syy(FrakFlag)=pressure(FrakFlag);
Sxy(FrakFlag)=0;

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
 

Option='C'; 

bottomboundflag=MidPoint(:,2)==-H;
Fdisp=bottomboundflag;
Option='H'; %No displacement along bottom boundary (Fdisp flag) 

figure;
scatter(MidPoint(:,1),MidPoint(:,2),15,Sxx,'filled');


smallturnNo=4; %Number of els we don't let turn significantly at start.
if i<=smallturnNo
    WidenSweep=0;
end
if WidenSweep==1
    ThetaIn=linspace(-20,20,7);  %The list of directions we are going to check (degrees from the current tip angle)
elseif WidenSweep==2
    ThetaIn=linspace(-90,90,21); 
else
    WidenSweep=0;
    ThetaIn=linspace(-7,7,8);  %Not so much needed once the first orientation is defined (defaults to larger searchs if this is not wide enough)
end
noNewels=1;
%% Find end loc.
[EndElsLoc] = GetCrackTipElements2d(P1,P2);
EndElsFlg=EndElsLoc>0;
if Walls==1
    %Find fracture tip - Edge point but lower that fracture walls
    Val=EndElsLoc./EndElsLoc; %ones for end els, nans for others.
    [~,I]=min(Val.*MidPoint(:,2));
    if MidPoint(I,2)==-H+HalfLength(I) %discount base of fracture
        Val(I)=NaN;
    end
    [~,I]=min(Val.*MidPoint(:,2));
else
    % %Get highest point
    [~,I]=max(MidPoint(:,2));
end
Indx=I;    
FuncArray=[H,g,rhoi,rhow,f];

figure;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(Option,'A')
    [Ds,Dn,Sxx,Syy,Sxy,StrainEnergy,InfMat]=SlipCalculator2dXtraCond(MidPoint,HalfLength,Sxx,Syy,Sxy...
        ,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option,...
        P1,P2,StrainEnergyDir,ThetaIn,Indx,noNewels,FuncArray);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Plot displacement discontinuity magnitude for all elements
GraphSlipVsElNum( Dn,Ds )

% Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Animate Fault movement.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %To do this look at file Animate2d.m


[K1,K2] = StressIntensity2d(EndElsLoc,Dn,Ds,HalfLength,E,nu,LineNormalVector,MidPoint,P1,P2);

%% Finding the current tip of the fracture
if Walls==1
    %Find fracture tip - Edge point but lower that fracture walls
    Val=EndElsLoc./EndElsLoc; %ones for end els, nans for others.
    [~,I]=min(Val.*MidPoint(:,2));
    if MidPoint(I,2)==-H+HalfLength(I) %discount base of fracture
        Val(I)=NaN;
    end
    [~,I]=min(Val.*MidPoint(:,2));
else
    % %Get highest point
    [~,I]=max(MidPoint(:,2));
end
Indx=I;    
Ktotal=sqrt(K1(I)^2+abs(K2(I))^2);

if EndElsLoc(I)==1
    EndPnt=[P1(I,1),P1(I,2)];
    FeM2Ev=normr([(P1(I,1)-MidPoint(I,1)),(P1(I,2))-MidPoint(I,2)]);
else %must be 2      

    EndPnt=[P2(I,1),P2(I,2)];
    FeM2Ev=normr([P2(I,1)-(MidPoint(I,1)),P2(I,2)-(MidPoint(I,2))]);
end

if K1(I)<0
    disp('fracture intepentration of the leading tip, leaving loop')
    i=i-1; %For plotting as we haven't added new el
    break
end

%% Computing K at the tip
ItpEls=sum(Dn(FrakFlag)<0);
if ItpEls>0
    fprintf(2,'\nFracture has X interpenetrating elements but continuing with prop loop: don''t trust results. X=\n')
    disp(ItpEls)
    pause(5)
end

%Compute with K ratio
if StrainEnergyDir==0
    %% If you want the prop dir from K
    K1(I)
    K2(I)
    K1(I)/(rhoi*g*H^(3/2))
    K2(I)/(rhoi*g*H^(3/2))
    Kc/(rhoi*g*H^(3/2))
    
    figure;
    scatter(K1(I),lenZH)
    figure;
    PlotFracture( P1(~EndElsFlg,:),P2(~EndElsFlg,:),'r' )
    hold on
    %Draw end els
    PlotFracture( P1(EndElsFlg,:),P2(EndElsFlg,:),'b' )
    hold off    
    
    %thetaNw=deg2rad(45);
    [thetaNw]=FindPropAngle(K1(I),K2(I));

    Ktotal=sqrt(K1(I)^2+abs(K2(I))^2);
 
else %Use strain energy
    figure; hold on
    %% If you are going to pick the best new prop direction using G (1 tip):
    StrainEnergyNw=-StrainEnergy;%fix sign
    scatter((ThetaIn),StrainEnergyNw,20,'k','filled');
    [mnstrEng,Indx_2]=min(StrainEnergyNw);
    thetaNw=deg2rad(ThetaIn(Indx_2));   
    scatter(ThetaIn(Indx_2),StrainEnergyNw(Indx_2),20,'b','filled')
    KtotalStrE=sqrt((abs(mnstrEng)*(E/(1-nu^2))));
    KtotalStrE2=sqrt((abs(mnstrEng)*(2*mu))/(1-nu));

    %Interpolate for smoothness?
    ThtHgh=linspace(min(ThetaIn),max(ThetaIn),500);
    vq2 = interp1((ThetaIn),StrainEnergyNw,ThtHgh,'spline');
    [mnstrEngItp,IndxItp]=min(vq2);
    thetaNw=deg2rad(ThtHgh(IndxItp));
    scatter(ThtHgh,vq2,5)
    scatter(ThtHgh(IndxItp),vq2(IndxItp),20,'r','filled')
    Ktotal=sqrt((abs(mnstrEngItp)*(E/(1-nu^2))));
    %Using str int approx
    Ktotal=sqrt(K1(I)^2+abs(K2(I))^2);

    if Indx_2==numel(StrainEnergy) || Indx_2==1
        
        if WidenSweep==2
            disp('sweep angle still not enuff, widensweep=-90:+90')
            disp('leaving loop because I can''t find a decent direction to grow. i=')
            disp(i)
            disp('K/K_c=')
            disp(Ktotal/Kc)

            %Checking inf mat stability... will error if not good
            try 
                InfMatrixCheck(InfMat); 
            catch
                disp('Your matrix was unstable, thats why...')
            end

            [Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d...
            (MidPoint(:,1),MidPoint(:,2),MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);
            figure;
            scatter(MidPoint(:,1),MidPoint(:,2),15,Uy)
            title('Vertical displacement')

            disp(K1(I))
            if K1(I)<0
                disp('The leading tip is also intepenetrating!!')
            end
            %[thetaNw]=FindPropAngle(K1(I),K2(I));
            i=i-1;
            break
        elseif WidenSweep==1
            disp('SWEEP ANGLE STILL NOT ENOUGH -> leaving loop')
            WidenSweep=2;
            i=i-1;
            continue            
            %break
        elseif WidenSweep==0
            if i>smallturnNo
                %restart loop with bigger search angle 
                WidenSweep=1;
                i=i-1;
                continue
            else
                WidenSweep=0;
            end
        end

    else
        WidenSweep=0; %Small search range
    end


end


%% Check we exceed Kc
Kc=0.1e6;

if Ktotal<Kc
    disp('leaving loop, K<Kc. i=')
    disp(i)
    disp(' K<Kc=')
    disp(Ktotal<Kc)

    i=i-1; %For plotting as we haven't added new el
    %disp('Break off!')
    break
end
KtotalLp(i)=Ktotal;
EndPointLp(i,:)=EndPnt;
KILp(i)=K1(I);
K2Lp(i)=K2(I);
if StrainEnergyDir==1
    GLp(i)=mnstrEng;
end

%% 3. Add new element with length & dir onto el list
StepSz=HalfLength(I)*2; %mean(HalfLength)*2; 
[FeM2Ev(1),FeM2Ev(2)] = RotateObject2d(FeM2Ev(1),FeM2Ev(2),thetaNw);
NewEnd=EndPnt+(FeM2Ev*StepSz);

if NewEnd(2)>0
    close all
    break %leave loop before half-space error
end

MrFlag=0; 
if MrFlag>1
    XBEG=P1NW(:,1);
    YBEG=P1NW(:,2);
    XEND=P2NW(:,1);
    YEND=P2NW(:,2);
else
    XBEG=P1(:,1);
    YBEG=P1(:,2);
    XEND=P2(:,1);
    YEND=P2(:,2);
end

%First checking we should actually add this new point
if any(~isnan(NewEnd))
    %If the current endpoint sits in P1. We add a new el, put this
    %endpoint in P2 and the new point in P1. This should maintain the
    %normal direction of the fracture
    if EndElsLoc(Indx)==1   %Do to BEG EL
        XBEG=[XBEG;EndPnt(1)];
        YBEG=[YBEG;EndPnt(2)];              

        XEND=[XEND;NewEnd(1)];
        YEND=[YEND;NewEnd(2)];
    else                    %Do to End EL
        XBEG=[XBEG;EndPnt(1)];
        YBEG=[YBEG;EndPnt(2)];              

        XEND=[XEND;NewEnd(1)];
        YEND=[YEND;NewEnd(2)];
    end
    [ MidPoint,HalfLength,P1,P2,LineNormalVector ]...
        = MidPoint_Orientation( XBEG,XEND,YBEG,YEND );            
end   

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

close all

%% current loop
figure; hold on
% Drawing fracture
%Normal els as red.
PlotFracture( P1,P2,'r' )
%Fixed elements shown as blue
id=find(Fdisp); 
PlotFracture(P1(id,:),P2(id,:),'b' )
%Titles
title('Predicted path using BEM: (σ_{θθ})max '), xlabel('x'), ylabel('z');axis equal
%Drawing Normals
%quiver(MidPoint(:,1),MidPoint(:,2),LineNormalVector(:,1),LineNormalVector(:,2))
h=(spring);
colormap(h)
scatter(EndPointLp(1:i,1),EndPointLp(1:i,2),15,KtotalLp(1:i)./Kc,'filled');
colorTitleHandle=colorbar;
%title(colorTitleHandle,'K/K_c')
ylabel(colorTitleHandle, 'K/K_c')
WhiteFigure;
axis('equal')
xlim([80 120])
ylim([-H -H+W])
caxis([0 25])
xlabel('x'), ylabel('z');
title('Current loop')

end
%% 
figure; hold on
% Drawing fracture
%Normal els as red.
PlotFracture( P1,P2,'r' )
%Fixed elements shown as blue
id=find(Fdisp); 
PlotFracture(P1(id,:),P2(id,:),'b' )
%Titles
title('Predicted path using BEM: (σ_{θθ})max '), xlabel('x'), ylabel('z');axis equal
%Drawing Normals
%quiver(MidPoint(:,1),MidPoint(:,2),LineNormalVector(:,1),LineNormalVector(:,2))
h=(spring);
colormap(h)
scatter(EndPointLp(1:i,1),EndPointLp(1:i,2),15,KtotalLp(1:i)./Kc,'filled');
colorTitleHandle=colorbar;
%title(colorTitleHandle,'K/K_c')
ylabel(colorTitleHandle, 'K/K_c')

%Now compute stresses without fracture and path based on most compressive stress
%direction:
X = linspace(0,W*2,grdsmpl);  
Y = linspace(-H+H*lenZH,-H+H,grdsmpl);
[X,Y] = meshgrid(X,Y);  %Location of observation points to compute stresses at.
[dimx,dimy] = size(X);  
[X,Y,Sxx,Syy,Sxy,Ux,Uy]=MainFrame2DFracture_StressesOnlyHNVD(H,W,L,nu,mu,noels,rhow,rhoi,g,f,X,Y,halfspaceflag,deltatau,Walls);

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

% Sxx=Sxx;
% Syy=Syy;

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

%Reshaping stresses to grid dimensions
lineobj =streamline(X,Y,S2x,S2y,W,-H+W*lenZH);
x=lineobj.XData;
y=lineobj.YData;
x(y>-1)=NaN;
y(y>-1)=NaN;

LineLen=0.25;
%figure;hold on
quiver(X(:),Y(:),S2x(:),S2y(:),LineLen,'k.')
plot(x,y,'b')
WhiteFigure;
axis('equal')
xlim([W-0.2*H W+0.2*H])
ylim([-H -H+W])
caxis([0 25])
xlabel('x'), ylabel('z');
title('Predicted path using BEM stress model: \sigma_2 parallel')

save('NEW_Panel_1250els_GCriterion_NoTractionNewEl_WithWalls_ZhStart=0_03_Dt=0_3_H=1000_W=1000_L=10000_f=0_7.mat')
function [ang]=FindPropAngle(K1,K2)
% %Right lateral movement so look in neg theta domain (see PF P375)
% if K2>0.0 
%     thetao = linspace(-0.5*pi,0.0,9000);
%     %Eq 9.78 Pollard:
%     angRad=(K1*sin(thetao))+(K2*((3.0*cos(thetao))-1.0));
% 
% %Left lateral movement
% else
%     thetao = linspace(0.0,0.5*pi,9000);
%     %Eq 9.78 Pollard:
%     angRad=(K1*sin(thetao))+(K2*((3.0*cos(thetao))-1.0));
% 
% end
% [mn,idx]=min(abs(angRad));
% ang=thetao(idx);

%Using the analytical formula
ang=-2*atan((-K1+sqrt(K1^2+8*K2^2))/(4*K2));

end