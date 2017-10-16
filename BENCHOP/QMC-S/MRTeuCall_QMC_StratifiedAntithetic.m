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
% By Josef H??k, josef.hook@it.uu.se, 2015
% Monte Carlo simulation of the Merton Jump model
% from analytical formula S = S0exp( (r  -lambda*xi - sigma^2/2)T + sigma W_t )*cumprod^N (1+Q)
% where 1+Q is lognormally distributed and N is the number of jumps (poisson) occured
% between [0,T]. 
% 
%
function U= MRTeuCall_QMC_StratifiedAntithetic(Sin,K,T,r,sigma, lambda, gamma,delta)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   

%randn('state', 100)

%{
%Bypass random generation step
%
% Set parameters
%
N = 1e6  % Num particles

% %%%%%%%%%%%%%%%%%%%%%%%
% Generate scrambled object 
% %%%%%%%%%%%%%%%%%%%%%%%
fprintf('Initializing scrambled faure sequence.');
SF =insfaur(1,N, 30); %scramble 30 digits
fprintf('...done\n');
%
Vals = zeros(1,N);

%
for i = 1:N
[Vals(i), SF] = gosfaur(SF);
end
%}
%save('Vals1e6.mat', 'Vals', 'N') 


load('Vals1e6.mat')

xi  = exp(gamma+delta.^2/2)-1;
mu = (r - lambda*xi - sigma.^2/2)*T;
Sigma = sigma*sqrt(T);

%
%
% Loads precalculated poisson Jumps if this has been performed previously
%
FileName = ['prodQplus1_', num2str(gamma), '_', num2str(delta),'_',num2str(2*N),'.mat'];
tic
if( exist(FileName, 'file') )
    load(FileName)
%disp('Loading poisson jumps')
else
    % Calculate number of Jumps
%  Should be possible to do this offline
Nt = poissrnd(lambda*T, 2*N,1);
prodQplus1 = zeros(2*N,1);
disp('Generating poisson jumps')
for k = 1:2*N
prodQplus1(k) = prod(lognrnd(gamma, delta, Nt(k),1));
end
toc
save(FileName, 'prodQplus1')
end


phi = normcdf(log(K./(S0.*prodQplus1)),mu,Sigma);
U1 = [Vals'; 1-Vals']; % antithetic
% Do not stratify the subset of particles
% that have zero probability 
F = phi + (1 -phi).*U1;

% Redraw particles that had zero probability
% We sometimes hit the machine precision 
% for normcdf. It can return 1 which correspond to
% a point a infty.
% We redraw all particles that hit this point
% and do not include them in the stratification
F(F==1) = U1(F==1);  % Reset
phi(F==1) = 0;       % Remove particles from stratification
% 

Z = norminv(F, mu, Sigma);
S = S0.*exp(Z).*prodQplus1;

% max operator almost not necessary here since we have drawn the 
% random numbers conditionally on S>K (stratification)
% except the rare events treated in the above section.
U(si) = exp(-r*T).*mean(max( S-K,0).*(1-phi));




end
%Test case: tic, MRTeuCallUI_QMC_StratifiedAntithetic(90,100,1,0.03,0.15, 0.4, -0.5, sqrt(0.16)), toc
%Xact = 7.542526012669795
%Error = abs(U-Xact)./abs(Xact)

U=U';