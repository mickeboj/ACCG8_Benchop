 %Copyright (C) 2015 Victor Shcherbakov

    %This file is part of BENCHOP.
    %BENCHOP is free software: you can redistribute it and/or modify
    %it under the terms of the GNU General Public License as published by
    %the Free Software Foundation, either version 3 of the License, or
    %(at your option) any later version.

    %BENCHOP is distributed in the hope that it will be useful,
    %but WITHOUT ANY WARRANTY; without even the implied warranty of
    %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %GNU General Public License for more details.

    %You should have received a copy of the GNU General Public License
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.

function U = HSTeuCall_RBFPUM(S,K,T,r,V,kap,th,sig,rho);


% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S - spot price
% K - strike price
% T - time to maturity
% r - risk free interest rate
% V - volatility
% kap - mean reversion
% th - mean level
% sig - vol diffusion
% rho - correlation
%
%---------  output  ----------------
% U - option value
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014

% ep=9.6, Nx=300, Ny=8, M=34, npx=12, npy=1, x,y = linspace, ymax=0.1


payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 9.6; %shape par

Nx = 300; % space points
Ny = 8;%ceil(Nx/4);
M = 35; % time step
npx = 12; % number of partitions
npy = 1;%ceil(npx/4);

op = '0';


% Generate suitable nodes for the type of option and domain
y_min = 0;
y_max = 0.1;
x_min = 0;
x_max = 4*K;
 
x = linspace(x_min,x_max,Nx);
y = linspace(y_min,y_max,Ny);
[xx,yy] = meshgrid(x,y);

xtype = zeros(size(xx)); % Interior
xtype(:,1) = 1; % X=0 boundary
xtype(:,end) = 2; % X=x_max boundary
xtype(end,2:end-1) = 3; % Y = y_max
xtype = xtype(:);

xc = [xx(:) yy(:)];

U = HestonRBFPU(S,K,T,r,V,kap,th,sig,rho,payoff,phi,ep,xc,M,npx,npy,op,xtype);
U = U';


%Ref values: 
% Uref = [2.302535842814927, 7.379832496149447, 14.974005277144057];

% Relative error:
% relerr = abs(U-Uref)./Uref

