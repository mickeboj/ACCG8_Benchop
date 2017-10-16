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

function [u]=Exact1D(x,t,sigma,r)
%
% The exact solution for the scaled equation with K=1.
% Seems correct when compared with the bc. (x-exp(-rt))
%  
pos = find(x==0);

x(pos)=x(pos)+1;

d1 = (log(x)+(r+0.5*sigma^2)*t)/(sigma*sqrt(t));
d2 = d1 - sigma*sqrt(t);

%
% From old file, it seems that N(x)=0.5*(1+erf(x/sqrt(2)))). Is this right?
%
Nd1 = 0.5*(1+erf(d1/sqrt(2)));
Nd2 = 0.5*(1+erf(d2/sqrt(2)));

u = x.*Nd1-exp(-r*t)*Nd2;
u(pos)=0;
