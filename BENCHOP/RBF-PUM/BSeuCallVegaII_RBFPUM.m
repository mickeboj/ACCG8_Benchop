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

function [Vega] = BSeuCallVegaII_RBFPUM(S,K,T,r,sig);

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
% Vega - Vega Greek
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014

% N = 5015, M = 1010, ep = 6, x = GetNode1D, ell = 0.25, np = N/80

S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 6; %shape par
N = 5015; % space points
M = 1010; % time step
np = ceil(N/80);


%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;
op = '0';

xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);


Vega=EuVegaRBFPU1D(S,K,T,r,sig,payoff,phi,ep,xc,M,np,op);
Vega = Vega';

%Ref values: 
% Vref=[10.689829936518105  12.307387008891816   0.224407208464969];

% Relative error:
% relerr = abs(Vega-Vref)./Vref



