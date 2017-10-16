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

function Vega = BSeuCallVegaI_RBFPUM(S,K,T,r,sig);

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

% N = 350, M = 160, ep = 0.1065, x = linspace, np = N/20


S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 0.1065; %shape par
N = 350; % space points
M = 160; % time step
np = ceil(N/20); % number of partitions
op = '0';


% Generate suitable nodes for the type of option and domain
x_min = 0;
x_max = 4*K;

xc = linspace(x_min,x_max,N)';

Vega=EuVegaRBFPU1D(S,K,T,r,sig,payoff,phi,ep,xc,M,np,op);
Vega = Vega';


%Ref values: 
% Vref=[32.770682446548165  38.413891530570481  28.995094522875153]; 

% Relative error:
% relerr = abs(Vega-Vref)./Vref

