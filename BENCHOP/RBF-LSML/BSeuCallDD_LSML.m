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

function U = BSeuCallDD_LSML(S,K,T,r,sig,D,alpha);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(S*(1-D)-K*exp(-r*t),0);  
phi = 'mq';
epf = 8;
epc = 12;
Nls = 600;
Nf = 0; % Only one level
Nc = 100;
M = 120;
op = '0';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
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
U=EuropeanLSML1D(S,K,T,r,sig,payoff,phi,epf,epc,xf,xc,xls,M,op)';
%
%Ref values: 
%Uref=[0.623811094545539 3.422712192881193 9.607009109943803]; 
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%end
