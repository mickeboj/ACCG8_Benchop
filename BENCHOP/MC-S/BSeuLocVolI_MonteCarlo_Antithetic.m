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
% European call with local vol model: dS = rSdt + sigma(S,t)SdW
% with MC 
% By Josef H??k, 2014
% Example usage: BSeuLocVolI_MonteCarlo_Antithetic(S0,K,T,r,@sigI)
%
%
function U = BSeuLocVolI_MonteCarlo_Antithetic(Sin,K,T,r,sigma)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   

% Set parameters
%
N = 1e5;  % Num particles
M = 1e4;  % Num time steps
t = linspace(0, T, M); 
dt = t(2)-t(1);

%sigma = @(s,t) 0.15 + 0.15*(0.5 + 2*t)*(s./100 - 1.2).^2./( (s./100).^2 + 1.44 );

% Rerun a larger simulation with N-N0 particles
S = ones(N, 1)*S0;
 
for l=1:M   
    % Antithetic pairs
    Z1 = randn(N/2,1 );
    Z2 = -Z1; 
    Z = [Z1; Z2];
    
    S = (1 + r*dt).*S + sigma(S,t(l)).*S.*sqrt(dt).*Z;
end

U(si) = mean(max(S-K, 0));
end
U=U';