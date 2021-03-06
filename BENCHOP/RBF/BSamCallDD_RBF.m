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

function U = BSamCallDD_RBF(S,K,T,r,sig,D,alpha);
% 
S = S(:);  
D = 0.03;
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
ep=28;
N=300;
M=60;
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
%
% Note that the double clustering does have an effect. 
%
xc = GetNodes1D(x_min,x_max,ntype,N,[K,(1-D)*K],ell,mode,scale);
%
U=AmericanDDRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep,alpha,D)';

%Ref values:
%Uref=[0.837358764004664 4.484034330080377 11.877216586990031]';
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%end
