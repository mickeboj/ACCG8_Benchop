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

function [U,relerr] = BSeuCallUII_RBF(S,K,T,r,sig);
% 
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
%

ep = 580;
N=4000; 
M = 1000;
op = '0';
varep = 'yes';
%
% Generate suitable nodes for the type of option and domain
ntype='cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;
xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);
%
%
U=EuropeanRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep)';
%
%Uref = 100*Exact1D(S/100,T,sig,r)
%Uref=[1.022192916849533e-04 0.033913177006145 0.512978189232605]';
% Relative error:
%relerr = (U-Uref)./Uref
