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

function U = BSupoutCallII_RBFPUM(S,K,T,r,sig,B);

% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S - spot price
% K - strike price
% T - time to maturity
% r - risk free interest rate
% sig - volatility
% B - barrier
%
%---------  output  ----------------
% U - option value
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014


S = S(:);  
payoff = @barrierpayoff;  
phi = 'mq';
ep = 5; %shape par
N = 3000; % space points
M = 1000; % time step
np = N/100; % number of partitions

%
% Generate suitable nodes for the type of option and domain
%
op = '0';
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = B; %1.25*K;
scale = K;

xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);

U=EuRBFPU1D(S,K,T,r,sig,payoff,phi,ep,xc,M,np,op);
U = U';


%Ref values: 
% Uref=[0.033913177006134   0.512978189232598   1.469203342553328]; 

% Relative error:
% relerr = abs(U-Uref)./Uref


function f = barrierpayoff(S,K,r,t)
  f = max(S-K*exp(-r*t),0);  
  % The barrier condition
  f(end)=0;


