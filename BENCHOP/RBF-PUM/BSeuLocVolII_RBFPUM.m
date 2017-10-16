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

function U = BSeuLocVolII_RBFPUM(S,K,T,r,sig);

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

% N = 130, M = 84, ep = 0.09, x = GetNodes1D, np = N/20

% N = 300, M = 150, ep = 1, x = GetNodes1D, np = N/16

% N = 150, M = 34, ep = 0.4, x = GetNodes1D, np = N/16  <--- best!

S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep = 0.4; %shape par
N = 150; % space points
M = 34; % time step
np = ceil(N/16); % number of partitions

%
% Generate suitable nodes for the type of option and domain
%
op = '0';
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;

xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);

for i = 1:length(S)
    U(i) = EuLocVolRBFPU1D(S(i),K,T,r,sig,payoff,phi,ep,xc,M,np,op);
end





