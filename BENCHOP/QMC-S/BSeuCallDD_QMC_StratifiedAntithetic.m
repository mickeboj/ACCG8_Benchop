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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.%
% By Josef H??k, josef.hook@it.uu.se
%
function U= BSeuCallDD_QMC_StratifiedAntithetic(Sin,K,T,r,sigma, D, alpha)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   
 

load('Vals1e6_3.mat')


% Run to dividend date
Ta = alpha*T;

mu1 = (r - sigma.^2/2)*Ta;
Sigma1 = sigma*sqrt(Ta);

Ua = [Vals2(:,1); 1-Vals2(:,1)]; % antithetic
%
% Generate GBM
%

Z = norminv(Ua, mu1, Sigma1);
S = S0.*exp(Z);


% Dividend
S = S - D*S;

% 
% Continue post dividend
mu = (r - sigma.^2/2)*(T-Ta);
Sigma = sigma*sqrt(T-Ta);
phi = normcdf( log(K./S),mu,Sigma);

U1 = [Vals2(:,2); 1-Vals2(:,2)];  % antithetic

% Do not stratify the subset of particles
% that have zero probability 
F = phi + (1 -phi).*U1;
% Redraw particles that had zero probability
% We sometime hit the machine precission 
% for normcdf. It can return 1 which correspond to
% a point at \infty.
% We redraw all particles that hit this point
% and do not include them in the stratification
F(F==1) = U1(F==1);  % Reset
phi(F==1) = 0;       % Remove particles from stratification
% 
Z = norminv(F, mu, Sigma);
S = S.*exp(Z);

U(si) = exp(-r*T).*mean(max(S-K,0).*(1-phi));
end
%Xact = 0.623811094545539
%Error = abs(U-Xact)./abs(Xact)
U=U';