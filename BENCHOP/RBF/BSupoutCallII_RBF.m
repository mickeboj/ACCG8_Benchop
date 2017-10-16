 %Copyright (C) 2015 Elisabeth Larsson

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

function U = BSupoutCallII_RBF(S,K,T,r,sig,B);
% 
S = S(:);  
payoff = @barrierpayoff;
phi = 'mq';
ep = 780;
N = 3000;
M = 2000;
op = '0';
varep = 'no';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = B;
scale = K;
xc = GetNodes1D(x_min,x_max,ntype,N,[K x_max],ell,mode,scale);
%
U=EuropeanRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep)';
%

%Ref values: 
%Uref=[0.033913177006134 0.512978189232598 1.469203342553328]'; 
% Relative error:
%relerr = (U-Uref)./Uref

function f = barrierpayoff(S,K,r,t)
  f = max(S-K*exp(-r*t),0);  
  % The barrier condition
  f(end)=0;