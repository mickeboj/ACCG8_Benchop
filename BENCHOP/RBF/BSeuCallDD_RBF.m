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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.% This gives the right answer.
function U = BSeuCallDD_RBF(S,K,T,r,sig,D,alpha);
% 
% The fastest way to solve the European call with parameter set 1 is to 
% cluster nodes around the strike and use the same shape of RBF everywhere.
%
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(S*(1-D)-K*exp(-r*t),0);  
phi = 'mq';
ep = 44;
N = 484;
M = 100;
op = '0';
varep = 'no';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;
xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);
%
U=EuropeanRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep)';
%
%Ref values:
%Uref=[0.623811094545539 3.422712192881193 9.607009109943803]';
% Relative error:
%relerr = (U-Uref)./Uref