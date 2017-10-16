 %Copyright (C) 2015 Erik Lindström

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

function [pExact,pDelta,pGamma,pVega]=BSCall_exact(S,K0,T,r,par)

% assume par is sigma, not sigma^2
sig=par;
sig2=sig*sig;

d1=(log(S/K0)+(r+.5*sig2)*T)/(sig*sqrt(T));
d2=d1-sig*sqrt(T);

pExact=S.*normcdf(d1)-K0*exp(-r*T)*normcdf(d2);
pDelta=normcdf(d1);
pGamma=normpdf(d1)./(S*sig*sqrt(T));
pVega=S.*normpdf(d1)*sqrt(T);


