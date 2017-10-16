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

function [U] = BSeuCallspread_RBF(S,T,r,sig1,sig2,rho);
% S is a 2xN vector  
if (size(S,2)>2)
  S = S';  
end
K = 0;
payoff = @(S1,S2,K,r,t) max(S1-S2-K*exp(-r*t),0);  
phi = 'mq';

ep = 40.75; 
N1 = 180;
N2 = 50; 
M=40; 

op = '0';
varep='no';
%
s2max = 150;
s1max = 4*(K+s2max);
[xc,xtype] = GetNodesSpread(s2max,N2,s1max,N1,K);
%
U=EuropeanRBF2D(S,K,T,r,sig1,sig2,rho,payoff,phi,ep,xc,xtype,M,op,varep)';

%Ref values: 
%Uref=[12.021727425647768, 5.978528810578943, 2.500244806693065, ...
%      2.021727425647768, 12.500244806693061]; 
%Uref=ExactSpread(T,S(:,1),S(:,2),sig1,sig2,rho);
% Relative error:
%  relerr = (U(:)-Uref(:))./Uref(:)
