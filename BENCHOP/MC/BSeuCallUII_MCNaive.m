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

function U= BSeuCallUII_MCNaive(Sin,K,T,r,sigma)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   
%randn('state', 100)


%
% Set parameters
%
N = 8e7;  % Num particles
M = 1024;

dt = T/M;


t = 0; 
S = S0;
while t<T
    S= S+ r*S*dt + sigma*S*sqrt(dt).*randn(N,1);
    t = t+dt;
end

U(si) = exp(-r*T).*mean(max(S-K,0));



end

%Xact = 2.758443856146076
%Error = abs(U-Xact)./abs(Xact)
U=U';