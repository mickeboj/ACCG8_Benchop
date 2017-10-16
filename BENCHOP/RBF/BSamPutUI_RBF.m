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

function U = BSamPutUI_RBF(S,K,T,r,sig);
% The code expects a column vector, so if S is in a row, fix this
S = S(:);  
payoff = @(S,K,r,t) max(K-S,0);  
phi = 'mq';
ep = 68;
N = 400;
M = 10000;
op = '0';
varep = 'no';
epen = 8e-7;
%
% Generate suitable nodes for the type of option and domain
%
ntype = 'cluster';
ell = 0.25;
mode = 'sinh';
x_min = 0;
x_max = 4*K;
scale = K;
xc = GetNodes1D(x_min,x_max,ntype,N,[0.82*K K],ell,mode,scale);
%
U=AmericanRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep,epen)';
%
% Relative error:

%Uref = [10.726486710094511
%     4.820608184813253
%     1.828207584020458];

%if (length(U)==length(Uref))
%  relerr = (U-Uref)./Uref
%end