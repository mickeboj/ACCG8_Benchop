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

function Delta = BSeuCallDeltaII_RBFPUM(S,K,T,r,sig);
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
% Delta - Delta Greek
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014


S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 4;
N = 3800; % space point
M = 300; % time step
varep = 'no';
np = N/200; % number of partitions
op = '1';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode='sinh';
x_min=0;
x_max = 4*K;
scale = K;
xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);

Delta=EuRBFPU1D(S,K,T,r,sig,payoff,phi,ep,xc,M,np,op);
Delta = Delta';

%Ref values: 
% Dref=[0.138001659888508   0.831964783803436   0.998616182178259];
% Relative error:
% relerr = abs(Delta-Dref)./Dref




