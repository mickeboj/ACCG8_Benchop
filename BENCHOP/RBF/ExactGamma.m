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

function [u]=ExactGamma(S,K,t,r,sig)
%
% The exact solution for the scaled equation with K=1.
% Seems correct when compared with the bc. (x-exp(-rt))
%  
pos = find(S==0);

S(pos)=S(pos)+1;

d1 = (log(S/K)+(r+0.5*sig^2)*t)/(sig*sqrt(t));
d2 = d1 - sig*sqrt(t);

%
% From old file, it seems that N(x)=0.5*(1+erf(x/sqrt(2)))). Is this right?
%
phid1 = 1/sqrt(2*pi)*exp(-0.5*d1.^2);
phid2 = 1/sqrt(2*pi)*exp(-0.5*d2.^2);

u = phid1./S/sig/sqrt(t);
u(pos)=0;
