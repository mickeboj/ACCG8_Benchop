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

function U= BSeuCallspread_MonteCarlo_Antithetic(Sin,T,r,sig1, sig2, rho)
Sin=Sin';
Nsin = size(Sin,2);
U =zeros(Nsin,1); 
for si = 1:Nsin

 S0_1 = Sin(1,si);   
 S0_2 = Sin(2,si);   
%randn('state', 100)

%S1-S2-K

%
% Set parameters
%
N = 1e8;  % Num particles

K = 0;

mu_1 = (r - sig1.^2/2)*T;
Sigma_1 = sig1*sqrt(T);
mu_2 = (r - sig2.^2/2)*T;
Sigma_2 = sig2*sqrt(T);


Z1t = randn(N/2,1);
Z1 = [Z1t; -Z1t];     %Antithetic
S1 = S0_1.*exp(mu_1 + Sigma_1.*Z1);
Z3t = randn(N/2,1);
Z3 = [Z3t; -Z3t]; % antithetic
Z2 = rho.*Z1 + sqrt(1-rho.^2).*Z3;
S2 = S0_2.*exp( mu_2 + Sigma_2.*Z2);

U(si) = exp(-r*T).*mean( max(S1-S2-K,0)); 


end

%Xact =5.978528810578943%12.021727425647768
%Error = abs(U-Xact)./abs(Xact)
U=U';