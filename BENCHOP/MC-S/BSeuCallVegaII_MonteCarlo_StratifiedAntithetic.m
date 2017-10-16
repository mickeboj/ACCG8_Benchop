 %Copyright (C) 2015 Josef Höök

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

function Vega= BSeuCallVegaII_MonteCarlo_StratifiedAntithetic(Sin,K,T,r,sigma)
Nsin = length(Sin);
Vega =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   

%
% Set parameters
%
N = 5e7;  % Num particles
% Antithetic pairs
%Z1 = randn(1,N/2 );
%Z2 = -Z1; 
%Ztot = [Z1, Z2];
%S = S0.*exp( (r - sigma.^2/2)*T+sigma*sqrt(T)*Ztot );

mu = (r - sigma.^2/2)*T;
Sigma = sigma*sqrt(T);
phi = normcdf(log(K/S0),mu,Sigma);
U0 = rand(N/2,1);
U1 = [U0; 1-U0]; % antithetic
F = phi + (1 -phi).*U1;
Z = norminv(F, mu, Sigma);
S = S0.*exp(Z);

Vega(si) = exp(-r*T).*mean( ( (log(S./S0) - (r + sigma^2/2)*T)./sigma ).*S )*(1-phi) ;
end

%test:
%tic,BSeuCallVegaII_MonteCarlo_StratifiedAntithetic(97,100,0.25,0.1,0.01),toc
%Xact = 10.689829936518105 
%Error = abs(Vega-Xact)./abs(Xact)

Vega=Vega';

