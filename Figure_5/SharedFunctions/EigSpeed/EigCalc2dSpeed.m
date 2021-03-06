function [S1,S2,S1dir,S2dir]=EigCalc2dSpeed(Sxx,Syy,Sxy)
% EigCalc2dSpeed: Calculates the 2D principal stress/strain magnitudes and 
%                   directions from input tensors. Like the function 
%                   EigCalc2d but around twice as fast. Can be used for any
%                   symmetrical 2*2 tensor Matrix and will return results
%                   of directions and magnitudes faster than just calling
%                   Eig. 
%   
% usage #1: For stress:
% [S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy)
%
% usage #2: For strain:
% [E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);
%
% Arguments: (input)
% Sxx,Syy,Sxy 		- The stress tensor components at a point, each can be a
%                    column vector.
%
% Arguments: (output)
% S1,S2       		- Principal stress component magnitudes Sigma 1 and
%                    Sigma 2.
%
% S1dir,S2dir       - Principal stress directions (direction cosines). Each
%                    will be a n*2 column vector [CosAx,CosAy] of this
%                    direction. 
%
% Example usage 1:
%
% %Calculating directions for a 2D stress tensor
% Sxx=0.2; Syy=-1.5; Sxy=1; 
% [S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy)
% %Drawing
% DrawS1S2Directions( 0,0,S1dir,S2dir )
%
% Example usage 2: Timing its speed compared to default calc
%
% n=1000000; %1000000;
% Sxx=rand(n,1);
% Syy=randn(n,1)*2;
% Sxy=ones(n,1);
% tic
% %Old eig func
% [S1,S2,S1dir,S2dir] = EigCalc2d(Sxx,Syy,Sxy);
% disp('Done part 1')
% toc
% tic
% %New fast func
% [E1,E2,E1dir,E2dir] = EigCalc2dSpeed(Sxx,Syy,Sxy);
% disp('Done part 2')
% toc
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Preallocating array
n = numel(Sxx);
tensor=zeros (2,2,n);

%Filling this array with 2x2 tensors, each 3rd dim is each point
%Accumulated using indexing
%[Sxx,Sxy]
%[Syx,Syy]
tensor(1,1,1:1:end) = Sxx(1:1:end,:);
tensor(1,2,1:1:end) = Sxy(1:1:end,:);
tensor(2,1,1:1:end) = Sxy(1:1:end,:);
tensor(2,2,1:1:end) = Syy(1:1:end,:);

%Eig can't handle nan's so we turn these to 0's and put the calculated s1s2 to nans after 
NanFlag = isnan(tensor);
tensor(NanFlag)=0;
NanFlag = isinf(tensor);
tensor(NanFlag)=0;

%Do the calculation
D = eig2(tensor); %see base of script and functions it calls
%Sort results
D = sortrows(D);
%Transpose
D=D';
%Extract
S2 = D(:,1);
S1 = D(:,2);

%Create vars for loop (preallocate)
S1dir=zeros(2,n);
%Create Identity matrix
Ident = repmat(eye(2),1,1,n);

for i=1:n
    
    %Equation we want to solve
    %[A-(Lambda * Identity)]v=0
    
    %Grab the current bit in the loop
    I=Ident(:,:,i);
    A=tensor(:,:,i);
    
    %[Lambda * Identity] in eigenvector equation
    One=S1(i)*I;
    
    %[A-(Lambda * Identity)]
    vect1=A-One;
    
    %[A-(Lambda * Identity)]v=0. Finding vector [v]
    %https://de.mathworks.com/help/matlab/ref/decomposition.html
    [Q1,~] = qr(vect1,0);

    %Extracting results
    S1dir(:,i)=(Q1(:,2));
    
end

S1dir=S1dir';

%We know S2 is perpendicular to S1. 
S2dir=[S1dir(2,:),-S1dir(1,:)];


function D = eig2(A)
% function D = eig2(A)
%
% Copyright (c) 2010, Bruno Luong
% All rights reserved.
%
%
% Compute in one shot the eigen-values of multiples (2 x 2) matrices
%
% INPUT:
%   A: (2 x 2 x n) array
% OUTPUT:
%   D: (2 x n). EIG2 returns in D(:,k) three eigen-values of A(:,:,k)
%
% See also: ParabolaRoots, eig3, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 27-May-2010

if size(A,1) ~= 2 || size(A,2) ~= 2
    error('A must be [3x3xn] array');
end

A = reshape(A, 4, []).';

P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,4));

% Determinant
P1 = A(:,1).*A(:,4) - A(:,2).*A(:,3);

% Find the roots of characteristic polynomials
D = ParabolaRoots(P3, P2, P1).';

end % eig2

end