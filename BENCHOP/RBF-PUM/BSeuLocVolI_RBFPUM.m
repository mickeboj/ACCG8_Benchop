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

function U = BSeuLocVolI_RBFPUM(S,K,T,r,sig);

% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S - spot price
% K - strike price
% T - time to maturity
% r - risk free interest rate
% sig - local volatility funtion
%
%---------  output  ----------------
% U - option value
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014

% ep = 0.1, x = linspace, N = 250, M = 30, np = N/20


S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 0.1; %shape par
N = 250; % space points
M = 30; % time step
np = ceil(N/20); % number of partitions
op = '0';

% Generate suitable nodes for the type of option and domain
x_min = 0;
x_max = 4*K;

xc = linspace(x_min,x_max,N)';

sig0 = @(r,S0,s,t) sig(s,t); 

U=EuLocVolRBFPU1D(S,K,T,r,sig0,payoff,phi,ep,xc,M,np,op);
U = U';

%Ref values: 
% Uref=[2.95517, 7.64217, 14.78425]; 

% Relative error:
% relerr = abs(U-Uref)./Uref

