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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.% I could compare three levels with BSMain, just to see if this becomes
% important after a while. This will be interesting to see.
function [U,err] = BSeuCallUII_LSML(S,K,T,r,sig);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';

epf = 44;
epc = 32;
Nls = 700;
Nf = 0;
Nc = 220;
M = 500;


op = '0';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
%ntype = 'uni';
ell = 0.25;
%mode = 'cheb';
mode = 'sinh';
x_min = 60;
x_max = 160;%4*K;
scale = K;
locs = [K x_max];
xc = GetNodes1D(x_min,x_max,ntype,Nc,locs,ell,mode,scale);
xf = GetNodes1D(x_min,x_max,ntype,Nf,locs,ell,mode,scale);
xls = GetNodes1D(x_min,x_max,ntype,Nls,locs,ell,mode,scale);
%
U=EuropeanLSML1D(S,K,T,r,sig,payoff,phi,epf,epc,xf,xc,xls,M,op)';
%
%Ref values: 
%Uref=[0.033913177006141 0.512978189232598 1.469203342553328]; 
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%  err=relerr;
%else
%  err = U-K*Exact1D(S(:)'/K,T,sig,r);
%end
