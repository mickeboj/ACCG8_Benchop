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

function [U] = BSeuCallUI_RBFPUM(S,K,T,r,sig);
% 
% The fastest way to solve the European call with parameter set 1 is to 
% cluster nodes around the strike and use the same shape of RBF everywhere.
%
% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S - spot price
% K - strike price
% T - time to maturity
% r - risk free interest rate
% sig - volatility
%
%---------  output  ----------------
% U - option value
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014


S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 0.2; %shape par
N = 180; % space points
M = 40; % time step
np = 8; % number of partitions
op = '0';

%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.2;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;

xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);

[U]=EuRBFPU1D(S,K,T,r,sig,payoff,phi,ep,xc,M,np,op);
U = U';

%Ref values: 
% Uref=[2.758443856146076   7.485087593912603   14.702019669720769]; 

% Relative error:
% relerr = abs(U-Uref)./Uref

