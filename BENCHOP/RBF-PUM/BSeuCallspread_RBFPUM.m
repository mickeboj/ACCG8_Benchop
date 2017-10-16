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

function U = BSeuCallspread_RBFPUM(S,T,r,sig1,sig2,rho);

% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S - spot price
% T - time to maturity
% r - risk free interest rate
% sig1 - volatility on 1st asset
% sig2 - volatility on 2nd asset
% rho - correlation
%
%---------  output  ----------------
% U - option value
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014

% ep=0.04, Nx=240, Ny=Nx/4, M=30, npx=16, npy=npx/4, x = linspace, 
% ymax=200, xmax = 5*ymax

% ep=0.04, Nx=129, Ny=Nx/5, M=24, npx=12, npy=npx/4, x = linspace, 
% ymax=150, xmax = 4*ymax    <---best!!!


payoff = @(S,K,r,t) max(S(:,1)-S(:,2),0);  
phi = 'mq';
ep = 0.04; %shape par
Nx = 129; % space points
Ny = ceil(Nx/5);
M = 30; % time step
npx = 12; % number of partitions
npy = ceil(npx/4);

op = '0';
K = 0;

% Generate suitable nodes for the type of option and domain
y_min = 0;
y_max = 150;
x_min = 0;
x_max = 4*y_max;
x = linspace(x_min,x_max,Nx);
y = linspace(y_min,y_max,Ny);
[xx,yy] = meshgrid(x,y);
xc = [xx(:) yy(:)];

U=EuRBFPU2Dspread(S,K,T,r,sig1,sig2,rho,payoff,phi,ep,xc,M,npx,npy,op);
U = U';




