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

function U = BSeuLocVolII_LSML(S,K,T,r,sig);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';

epf = 10;
epc = 8;
Nls = 420;
Nf = 0; % This should mean nothing happening at this level
Nc = 84;
M = 70;

op = '0';
limited='yes';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'cheb';
x_min = 0;
x_max = 4*K;
scale = K;
locs = [x_max];
xc = GetNodes1D(x_min,x_max,ntype,Nc,locs,ell,mode,scale);
xf = GetNodes1D(x_min,x_max,ntype,Nf,locs,ell,mode,scale);
xls = GetNodes1D(x_min,x_max,ntype,Nls,locs,ell,mode,scale);
%
for k=1:length(S)
  S0=S(k);
  sigS0 = @(s,t) sig(r,S0,s,t); 
  U(1,k)=EuropeanLocVolLSML1D(S0,K,T,r,sigS0,payoff,phi,epf,epc,xf,xc,xls,M, ...
                              op,limited);
end

%
%Ref values: 
%Uref=[2.1433 6.8740 14.2983]; 
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%end
