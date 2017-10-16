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
% Stratified sampling of an up and out Call 
% Reference: Conditionong on one-step survival for barrier option
% simulations
% 

function [U]=BSupoutCallII_MonteCarlo_StratifiedAntithetic(Sin,K,T,r,sigma,B)

Nsin = length(Sin);
U =zeros(Nsin,1); 

N = 7e7;  % Num particles
M = 1e2;  % Intermediate time steps
t = linspace(0,T, M);
dt = t(2);

for si = 1:Nsin

  S = Sin(si).*ones(N,1);
    mu = (r - sigma.^2/2)*dt;
    Sigma = sigma*sqrt(dt);

    L = ones(N,1); % Cumulated likelihood weights
 for l=1:M
 
     phi = normcdf(log(B./S),mu,Sigma);
U0 = rand(N/2,1);
U1 = [U0; 1-U0]; % antithetic
F = phi.*U1;
Z = norminv(F, mu, Sigma);
Sold = S;
S = S.*exp(Z);

% Check if the path has crossedbetween S_i-1 and S_i
 Pl = exp( -2*(B-Sold).*(B-S)./( sigma^2.*Sold.^2*dt ) );
 L(Pl>=rand(N,1)) = 0; % kill paths that     
     
%{    phi_old = normcdf(log(B./S),mu,Sigma);
% phi = normcdf( (log(B./S) + mu)./Sigma);
%    err_phi = norm(phi_old -phi)
%    U0 = rand(N/2,1);
%    U1 = [U0; 1-U0]; % antithetic
%    F = phi.*U1 + (1 -phi);
%    Z = norminv(F, mu, Sigma);
%    S = S.*exp(Z);   
  %}
    L = L.*(phi);
 end

%
%
Val = max(S-K,0)'*L;
U(si) = exp(-r*T).*Val/N;   

%Xact = 1.822512255945242 
%Error = abs(Xact-U(si))/abs(Xact)
end
U = U';

