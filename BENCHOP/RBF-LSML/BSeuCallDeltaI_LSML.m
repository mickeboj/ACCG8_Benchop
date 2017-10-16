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

function Delta = BSeuCallDeltaI_LSML(S,K,T,r,sig);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';
epf = 4;
epc = 0.5;
Nls = 380;
Nf = 40; 
Nc = 20;
M = 60;

epf = 4;
epc = 0.5;
Nls = 380;
Nf = 40; 
Nc = 20;
M = 60;
op = '1';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;
locs = [K x_max];
xc = GetNodes1D(x_min,x_max,ntype,Nc,locs,ell,mode,scale);
xf = GetNodes1D(x_min,x_max,ntype,Nf,locs,ell,mode,scale);
xls = GetNodes1D(x_min,x_max,ntype,Nls,locs,ell,mode,scale);
%
Delta=EuropeanLSML1D(S,K,T,r,sig,payoff,phi,epf,epc,xf,xc,xls,M,op)';
%
%Ref values: 
%Dref=[0.334542751969886 0.608341880846395 0.818694517094515];
% Relative error:
%if (length(Delta)==length(Dref))
%  relerr = (Delta-Dref)./Dref
%end
