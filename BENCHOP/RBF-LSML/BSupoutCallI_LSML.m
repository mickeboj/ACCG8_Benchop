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

function [U] = BSupoutCallI_LSML(S,K,T,r,sig,B);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @barrierpayoff;
phi = 'mq';
epf = 0;
epc = 4;
Nls = 280;
Nf = 0;%40; % This should mean nothing happening at this level
Nc = 32;
M = 80;
op = '0';
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'cheb';
x_min = 0;
x_max = B;
scale = K;
locs = [K x_max];
xc = GetNodes1D(x_min,x_max,ntype,Nc,locs,ell,mode,scale);
xf = GetNodes1D(x_min,x_max,ntype,Nf,locs,ell,mode,scale);
xls = GetNodes1D(x_min,x_max,ntype,Nls,locs,ell,mode,scale);
%
[U]=EuropeanLSML1D(S,K,T,r,sig,payoff,phi,epf,epc,xf,xc,xls,M,op)';
%
%Ref values: 
%Uref=[1.822512255945242   3.294086516281595   3.221591131246868]; 
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%end

function f = barrierpayoff(S,K,r,t)
  f = max(S-K*exp(-r*t),0);  
  % The barrier condition
  f(end)=0;