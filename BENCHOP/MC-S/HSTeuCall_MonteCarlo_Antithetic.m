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
%  Simulation of the Heston model using the QE scheme
%  Source: http://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/1826b88b152e65a7c12574b000347c74/$FILE/LeifAndersenHeston.pdf
%  [1]
% Test:  tic, HSTeuCall_MonteCarlo_Antithetic(90,100,1,0.03,0.0225,2,0.0225,0.25,-0.5), toc
%
function U=HSTeuCall_MonteCarlo_Antithetic(Sin,K,T,r,V0,kappa,theta,sigma,rho)
%function U=HSTeuCall_MonteCarlo_Antithetic()

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   



%
%   dS = rSdt +  sqrt(V)S dW_s
%   dV   =  kappa(theta-V)*dt + sigma*sqrt(V)dW_v
%

% parameters
dt = 1/256;

% Simulation specific
N = 6e7;  % num particles


V = V0*ones(1,N);%0.0225*ones(1,N);
S = S0*ones(1,N);
lnS = log(S);


% Algorithm params
PsiCrit = 1.5;
% For central difference 
gamma1= 1/2;
gamma2 = 1/2;

K0 = -rho.*kappa.*theta/sigma.*dt;
K1 = gamma1.*dt*( kappa.*rho/sigma -1/2 )  - rho/sigma;
K2 = gamma2.*dt*( kappa.*rho/sigma -1/2 )  + rho/sigma; 
K3 = gamma1.*dt*(1 - rho.^2); 
K4 = gamma2.*dt*(1 -rho.^2);

t = 0;
while t<T
    %
    % Step 1
    %
    m = theta + (V - theta).*exp(-kappa*dt);
    s2 = (V*sigma.^2.*exp(-kappa*dt))./kappa.*(1-exp(-kappa*dt)) + theta.*sigma^2./(2*kappa).*(1- exp(-kappa*dt)).^2;
    Psi = s2./m.^2;
    
    %
    % Step 2
    % Draw a uniform random number
    tmp = rand(1,N/2);
    Uv= [tmp, 1-tmp];
  
      % default mode
        % Use eq 26
        p = (Psi-1)./(Psi + 1);
        beta = (1-p)./m;
        Vdt = PsiInv(Uv, p, beta);
 
     % find cases where we use eq 23 instead
        idx = find(Psi<=PsiCrit);     
        b2 = 2./Psi(idx)-1 + sqrt(2./Psi(idx)).*sqrt(2./Psi(idx)-1);
        b = sqrt(b2);
        a = m./(1 + b2);
        
        Zv = norminv(Uv(idx));
        Vdt(idx) = a.*(b+Zv).^2;
        
    
    % 
    % S step
    % 
    tmp2 = randn(1,N/2);
    Z = [tmp2, -tmp2];
    lnSdt = lnS + K0 + K1*V + K2*Vdt + sqrt(K3*V + K4*Vdt).*Z + r*dt;
    
    
    V = Vdt;
    lnS = lnSdt;
    t = t+dt;
end

% Evaluate Payoff
U(si) = mean(exp(-r*T).*max(exp(lnS)-K,0));

% tic, HSTeuCall_MonteCarlo_Antithetic(90,100,1,0.03,0.0225,2,0.0225,0.25,-0.5), toc
%Xact = 2.302535842814927%7.379832496149447
%Error = abs(Xact- U(si))./abs(Xact)
end

U = U';

end

% Eq. 25 in [1]
function ret=PsiInv(U,p,beta)
% assumes 0<U<=1
ret = zeros(size(U));
ret(p<U) = log((1-p)./(1-U(p<U)))./beta;
end

