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

function U= BSeuCallDD_MonteCarlo_StratifiedAntithetic(Sin,K,T,r,sigma, D, alpha)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   

%randn('state', 100)
N = 9e7;

%
% Generate GBM
%

% Run to dividend date
Ta = alpha*T;

mu1 = (r - sigma.^2/2)*Ta;
Sigma1 = sigma*sqrt(Ta);
dW1 = randn(N/2,1);
dWa = [dW1; -dW1]; % antithetic
S = S0.*exp(mu1 + Sigma1.*dWa);

% Dividend
S = S - D*S;

% 
% Continue post dividend
mu = (r - sigma.^2/2)*(T-Ta);
Sigma = sigma*sqrt(T-Ta);
%dW2 = randn(N/2,1);
%dWb = [dW2; -dW2]; % antithetic
%Smin = min(S);
phi = normcdf( log(K./S),mu,Sigma);

U0 = rand(N/2,1);
U1 = [U0; 1-U0]; 
% Do not stratify the subset of particles
% that have zero probability 
F = phi + (1 -phi).*U1;

% Redraw particles that had zero probability
% We sometimes hit the machine precission 
% for normcdf. It can return 1 which correspond to
% a point a infty.
% We redraw all particles that hit this point
% and do not include them in the stratification
F(F==1) = U1(F==1);  % Reset
phi(F==1) = 0;       % Remove particles from stratification
% 
Z = norminv(F, mu, Sigma);
S = S.*exp(Z);





U(si) = exp(-r*T).*mean(max(S-K,0).*(1-phi));

end
%test: 
%tic, BSeuCallDD_MonteCarlo_StratifiedAntithetic(90,100,0.5,0.03,0.15,0.03, 0.8),toc
%Xact = 0.623811094545539
%Error = abs(U-Xact)./abs(Xact)
U=U';