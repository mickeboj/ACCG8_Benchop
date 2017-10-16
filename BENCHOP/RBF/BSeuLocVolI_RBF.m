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

function U = BSeuLocVolI_RBF(S,K,T,r,sig)
% Note that sig here is a function handle, not a constant. Send in @locvolI
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
% Dessa troligen tillr???ckliga. N???sta set skiljer i 6:e siffran.
ep = 14;
N = 154;
M = 30;
%ep = 48;
%N=600;
%M=200;
%ep=120;
%N=1000;
%M=1000;
%ep=100;
%N=1000;
%M=1000;
op = '0';
varep = 'no';
limited='no';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;
% This is not scale invariant... Hmm!
xc = GetNodes1D(x_min,x_max,ntype,N,K,ell,mode,scale);
%
U=EuropeanLocVolRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep,limited)';