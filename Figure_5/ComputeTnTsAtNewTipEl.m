function [Tn,Ts] = ComputeTnTsAtNewTipEl(MidPoint,LineNormalVector,FuncArray)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

H=FuncArray(1);
g=FuncArray(2);
rhoi=FuncArray(3);
rhow=FuncArray(4);
f=FuncArray(5);

Pw=f*(rhoi*g*H); %Water pressure at base of glacier

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
Pxx=pressure; %Pw
Pyy=pressure;
Pxy=0;

CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2);  
%Calculates the normal and shear tractions on the elements 
[ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,CosAx,CosAy );
[ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy);

end