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

function Vega= BSeuCallVegaII_QMC_StratifiedAntithetic(Sin,K,T,r,sigma)

Nsin = length(Sin);
Vega =zeros(Nsin,1);
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

%
% Generate GBM
%

mu = (r - sigma.^2/2)*T;
Sigma = sigma*sqrt(T);
phi = normcdf(log(K/S0),mu,Sigma);
U1 = [Vals'; 1-Vals']; % antithetic
F = phi + (1 -phi).*U1;

Z = norminv(F, mu, Sigma);
S = S0.*exp(Z);


Vega(si) = exp(-r*T).*mean( ( (log(S./S0) - (r + sigma^2/2)*T)./sigma ).*S )*(1-phi) ;

end
%Xact = 32.770682446548165
%Error = abs(Vega-Xact)./abs(Xact)

Vega=Vega';
