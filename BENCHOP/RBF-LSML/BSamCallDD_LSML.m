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

function U = BSamCallDD_LSML(S,K,T,r,sig,D,alpha);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';

epf = 4;
epc = 5;
Nls = 600;
Nf = 0; % Only using one level
Nc = 100;
M = 300;


op = '0';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'uni';
ell = 0.25;
mode = 'cheb';
x_min = 0;
x_max = 4*K;
scale = K;
locs = [K x_max];
xc = GetNodes1D(x_min,x_max,ntype,Nc,locs,ell,mode,scale);
xf = GetNodes1D(x_min,x_max,ntype,Nf,locs,ell,mode,scale);
xls = GetNodes1D(x_min,x_max,ntype,Nls,locs,ell,mode,scale);
%
[U] = AmericanDDLSML(S,K,T,r,sig,payoff,phi,epf,epc,xf,xc,xls,M,op,alpha,D)';
%
%Ref values: 
%Uref=[0.837358764004664 4.484034330080377 11.877216586990031]; 
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%end
