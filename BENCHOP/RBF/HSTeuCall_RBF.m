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

function U = HSTeuCall_RBF(S,K,T,r,V,kap,th,sig,rho);
S = S(:);  
payoff = @(S,K,r,t) max(S-K*exp(-r*t),0);  
phi = 'mq';

ep = 20;
Ns = 320;
Nv = 8;
M = 100;

op = '0';
%
% Generate uniform nodes in [0,smax]x[0,vmax]
%
smax = 4*K;
vmax = 0.1;

s = linspace(0,smax,Ns);
v = linspace(0,vmax,Nv);

[xx,yy] = meshgrid(s,v);
%
% Create information about which points are at which boundary
%
xtype = zeros(size(xx)); % Interior
xtype(:,1) = 1; % X=0 boundary
xtype(:,end) = 2; % X=x_max boundary

xtype = xtype(:);
xc = [xx(:),yy(:)];
%
U=HestonRBF2D(S,K,T,r,V,kap,th,sig,rho,payoff,phi,ep,xc,xtype,M,op)';
%

%Ref values: 
%Uref=[2.302535842814927, 7.379832496149447, 14.974005277144057];
% Relative error:
%if (length(U)==length(Uref))
%  relerr = (U(:)-Uref(:))./Uref(:)
%end
