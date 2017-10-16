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

function u=ExactSpread(t,S1,S2,sig1,sig2,rho)
%
% To have something to compare with.
%  
x = S1./S2;

pos = find(S1.*S2<=eps);
x(pos)=1;

sigma = sqrt(sig1^2 - 2*rho*sig1*sig2 + sig2^2);
  
d1 = (log(x)+(0.5*sigma^2)*t)/(sigma*sqrt(t));
d2 = d1 - sigma*sqrt(t);

%
% From old file, it seems that N(x)=0.5*(1+erf(x/sqrt(2)))). Is this right?
%
Nd1 = 0.5*(1+erf(d1/sqrt(2)));
Nd2 = 0.5*(1+erf(d2/sqrt(2)));

u = S1.*Nd1-S2.*Nd2;

u(pos)=0;
  