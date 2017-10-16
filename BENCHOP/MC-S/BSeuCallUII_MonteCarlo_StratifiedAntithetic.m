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

function U= BSeuCallUII_MonteCarlo_StratifiedAntithetic(Sin,K,T,r,sigma)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   
%randn('state', 100)


%
% Set parameters
%
N = 5e7;  % Num particles


mu = (r - sigma.^2/2)*T;
Sigma = sigma*sqrt(T);
phi = normcdf(log(K/S0),mu,Sigma);
U0 = rand(N/2,1);
U1 = [U0; 1-U0]; % antithetic
F = phi + (1 -phi).*U1;

Z = norminv(F, mu, Sigma);
S = S0.*exp(Z);
% max operator not necessary here since we have drawn the 
% random numbers conditionally on S>K (stratification)
U(si) = exp(-r*T).*mean(S-K)*(1-phi);



end
% test: tic, BSeuCallUII_MonteCarlo_StratifiedAntithetic(97,100, 0.25, 0.1,
% 0.01), toc
%Xact = 0.033913177006141
%Error = abs(U-Xact)./abs(Xact)
U=U';