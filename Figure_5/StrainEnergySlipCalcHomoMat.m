function [StrainEnergyOut]=StrainEnergySlipCalcHomoMat(D,StressInf,HalfLength,MidPoint,P1,P2,...
            halfspace,nu,E,Fdisp,TipIndxs,noNewels,ThetaIn,Tn,Ts,TnIn,TsIn,AddTractionToNewEls,FuncArray )

%% For the calculator
% nu/mu  % The elastic consts
% ----------Boundary conditions:
% Tn Ts %B= [Tn;Ts];  
% ----------How we test the new tips
% ThetaIn=[90,0];   %The list of directions we are going to check (degrees from the current tip angle)
% noNewels=25;      %The number of new tip elements (1 if using Macca style calc - matches K from BEM approx)
% doubledelflection=0;      %If we want the new tip to deflect in two directions
% AddTractionToNewEls=0;    %If we want the new elements traction free or -> not... (not for Macca style)
% Indx - tip index's you are interest in: Computed outside using:
[EndElsLoc] = GetCrackTipElements2d(P1,P2);
% EndElsFlg=EndElsLoc>0;
% [TipIndxs]=find(EndElsLoc~=0);
Indx=TipIndxs;

%Quick fix
TnIn=Tn;
TsIn=Ts;


%% STRAIN ENERGY RELEASE
[ Dn,Ds ] = ExtractArraysFromVector( D );
%TnNew at els with zero vert displacement boundary condition: TnNew(FD)
%should be ~ zero
TnNew =  -StressInf.DsTn*Ds + -StressInf.DnTn*Dn;	
TsNew =  -StressInf.DsTs*Ds + -StressInf.DnTs*Dn;	
StrainEnergy=sum(Dn.*TnNew.*HalfLength+Ds.*TsNew.*HalfLength);

StE=[];

%0 when straight -> supply vec if you want to test multiple...
Theta=ThetaIn;%[90,0];%linspace(-MxAng,MxAng,20);%MxAng*2+1);
Theta=deg2rad(Theta);
%noNewels=1; %Number of new tip elements  (free bound cond)
disp('Midway:NewTipSampling')
disp(noNewels)

%doubledelflection=0;


for J=1:numel(Theta)

    for i =1:numel(Indx)

        I=Indx(i); %Get current index    
        if EndElsLoc(Indx(i))==1
            EndPnt=[P1(I,1),P1(I,2)];
            FeM2Ev=normr([(P1(I,1)-MidPoint(I,1)),(P1(I,2))-MidPoint(I,2)]);
        else %must be 2      
            EndPnt=[P2(I,1),P2(I,2)];
            FeM2Ev=normr([P2(I,1)-(MidPoint(I,1)),P2(I,2)-(MidPoint(I,2))]);
        end

        thetaNw=Theta(J);

        %% 3. Add new element with length & dir onto el list
        StepSz=(min(HalfLength)*2)*noNewels; 
        %StepSz=StepSz/100; disp('way small!')
        %StepSz=(mean(HalfLength))*noels; 
        %noels=25; %Number of new tip elements  (free bound cond)
        disp('NewtipLength')
        newtiplength=(min(HalfLength)*2)*noNewels;
        disp(newtiplength)
        
        %If you want at current angle
        [FeM2Ev(1),FeM2Ev(2)] = RotateObject2d(FeM2Ev(1),FeM2Ev(2),thetaNw);
        NewEnd=EndPnt+(FeM2Ev*StepSz);

        XBEG=P1(:,1);
        YBEG=P1(:,2);
        XEND=P2(:,1);
        YEND=P2(:,2);

        %First checking we should actually add this new point
        if any(~isnan(NewEnd))
            %If the current endpoint sits in P1. We add a new el, put this
            %endpoint in P2 and the new point in P1. This should maintain the
            %normal direction of the fracture
            if EndElsLoc(Indx(i))==1 %Do to BEG EL
                %% new - for multiple els along new part...

                x=linspace(NewEnd(1),EndPnt(1),noNewels+1);
                y=linspace(NewEnd(2),EndPnt(2),noNewels+1);
                xbegnew=x(1:end-1)';
                xendnew=x(2:end)';
                ybegnew=y(1:end-1)';
                yendnew=y(2:end)';  
                XBEG=[xbegnew;XBEG];
                YBEG=[ybegnew;YBEG];
                XEND=[xendnew;XEND];
                YEND=[yendnew;YEND]; 
                FdispNW=[Fdisp;0]; %New el not locked
                addedontop=1;
                
            elseif EndElsLoc(Indx(i))==2 %Do to End EL
                %% 2.  new - for multiple els along new part...
                x=linspace(NewEnd(1),EndPnt(1),noNewels+1);
                y=linspace(NewEnd(2),EndPnt(2),noNewels+1);
                xbegnew=x(1:end-1);
                xendnew=x(2:end);
                ybegnew=y(1:end-1);
                yendnew=y(2:end);  
                XBEG=[XBEG;(xendnew )'];
                YBEG=[YBEG;(yendnew)'];
                XEND=[XEND;(xbegnew)'];
                YEND=[YEND;( ybegnew )'];    
                FdispNW=[0;Fdisp]; %New el not locked
                addedontop=0;
               
            end
            [ MidPointNW,HalfLengthNW,P1NW,P2NW,LineNormalVectorNW ]...
                = MidPoint_Orientation( XBEG,XEND,YBEG,YEND );            
        else
            disp(NewEnd)
            error('Why does NewEnd contain NAN vals?')
        end

    close 
    figure;
    PlotFracture( P1NW,P2NW,'r' ); axis('equal');
    %Drawing Normals
    quiver(MidPointNW(:,1),MidPointNW(:,2),LineNormalVectorNW(:,1),LineNormalVectorNW(:,2))
    %close

    %% Setup boundary condition
    %Add this extra el
    CosAx=LineNormalVectorNW(:,1);
    CosAy=LineNormalVectorNW(:,2);  
    %Making remote Sxx,Syy,Sxy a list the size of the number of elements
    ZeroList = zeros(numel(CosAx),1);
    %Zero tractio new el

    [TnTipEl,TsTipEl] = ComputeTnTsAtNewTipEl(MidPointNW(end,:),LineNormalVectorNW(end,:),FuncArray);

    if addedontop==1
        Tn= [TnTipEl;TnIn]; %Assuming the whole fracture has the same stress
        Ts= [TsTipEl;TsIn];
    else
        Tn= [TnIn;TnTipEl]; %Assuming the whole fracture has the same stress
        Ts= [TsIn;TsTipEl]; 
    end

    if AddTractionToNewEls==1
%         %Assuming linear grad
%         Tn= (MidPointNW(:,2)-b).*gamma; %Assuming the whole fracture has the same stress
%         Ts= ZeroList;  
    end

    %% Setup Inf mats
    [StressInfNw,DispInfNw]=CalculateInfluenceMatrices2d(halfspace,MidPointNW,HalfLengthNW,nu,E,LineNormalVectorNW,Fdisp );
    [StressInfNw,Tn,Ts,NUM]=FixingDisp_InfRowSwap2dHanwen(StressInfNw,DispInfNw,numel(FdispNW),Tn,Ts,FdispNW);


    Ats = [-StressInfNw.DnTn,-StressInfNw.DsTn];            %StressInf = rmfield(StressInf,{'DnTn','DsTn'});   
    Ash = [-StressInfNw.DnTs,-StressInfNw.DsTs];            
    %clear DsTs DnTs DsTn DnTn                                              
    ANW=  [Ats;Ash];
    %clear Ash Ats 


    %% Solve matrix system
    BNw= [ Tn;Ts ];
    DNw = ANW\BNw;
    [ DnNW,DsNW ] = ExtractArraysFromVector( DNw );
    
%     figure;hold on;axis('equal')
%     %Grab direction cosines of each element
%     CosAx=LineNormalVectorNW(:,1);
%     CosAy=LineNormalVectorNW(:,2);
%     %Compute the directions of the slip components.
%     DnCosine=LineNormalVectorNW;
%     DsCosine=[-CosAy,CosAx];
%     %Cartesian directions and magnitudes of slip
%     Dn_xy=(bsxfun(@times,DnNW,DnCosine));
%     Ds_xy=(bsxfun(@times,DsNW,DsCosine));
%     scatter(MidPointNW(:,1)+Dn_xy(:,1)+Ds_xy(:,1),...
%         MidPointNW(:,2)+Dn_xy(:,2)+Ds_xy(:,2));
%     scatter(MidPointNW(:,1)-Dn_xy(:,1)-Ds_xy(:,1),...
%         MidPointNW(:,2)-Dn_xy(:,2)-Ds_xy(:,2),'r');
%     scatter(MidPointNW(:,1),...
%         MidPointNW(:,2),'.k')
%     
%     xx=MidPointNW(end-noNewels:end,1)+Dn_xy(end-noNewels:end,1)+Ds_xy(end-noNewels:end,1);
%     yy=MidPointNW(end-noNewels:end,2)+Dn_xy(end-noNewels:end,2)+Ds_xy(end-noNewels:end,2);
%     maxx=max([abs(xx); abs(yy)]); %Max disp near tip
%     xlim([-maxx*2 maxx*2]);
%     ylim([-maxx*2 maxx*2]);
%     
%     hold off
    
    %Orig -no fixed disp els
    TnNew =  -StressInfNw.DsTn*DsNW + -StressInfNw.DnTn*DnNW;	
    TsNew =  -StressInfNw.DsTs*DsNW + -StressInfNw.DnTs*DnNW;	
    StrainEnergyNw=sum(DnNW.*TnNew.*HalfLengthNW+DsNW.*TsNew.*HalfLengthNW);
%     %Ignore els with displacement boundary condition: (~FD,:)
%     TnNew =  -StressInfNw.DsTn*DsNW + -StressInfNw.DnTn*DnNW;	
%     TsNew =  -StressInfNw.DsTs*DsNW + -StressInfNw.DnTs*DnNW;	
%     StrainEnergyNw=sum(DnNW(~FdispNw,:).*TnNew(~FdispNw,:).*HalfLengthNW(~FdispNw,:))+sum(DsNW.*TsNew.*HalfLengthNW);

    %Maccaferri 2010: 
    %Assuming we add just 1 el
    StrainEnergyChange=(StrainEnergyNw-StrainEnergy)/newtiplength; %*2 == see Dahm 2000 p.141 bottom left paragraph

    disp('theta')
    disp(Theta(J))
    disp('strain energy change')
    disp(StrainEnergyChange)  

    StE=[StE,StrainEnergyChange];   
    

    end

end

%Correct sign...
StrainEnergyOut=-StE;

end