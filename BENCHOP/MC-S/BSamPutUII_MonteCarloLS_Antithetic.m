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
% Implements Longstaff-Schwartz least squares method for pricing american options
% By Josef Hook
% Reference: http://www.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf
%
function U= BSamPutUII_MonteCarloLS_Antithetic(Sin,K,T,r,sigma)

Nsin = length(Sin);
U =zeros(Nsin,1);
for si = 1:Nsin

 S0 = Sin(si);   

%randn('state', 100)


%
% Set parameters
%
N = 5e1;  % Num particles
M = 1e3;  % Num time steps
t = linspace(0, T, M); 
dt = t(2)-t(1);

%
% Generate GBM
%
% Antithetic pairs
Z1 = randn(M,N/2 );
Z2 = -Z1; 
Ztot = [Z1, Z2];

St = exp( (r - sigma.^2/2)*dt+sigma*sqrt(dt)*Ztot );
S = cumprod([S0*ones(1,N); St],1); 
% Note that size will be of size M+1


%
% Code below have been validated 2014-12-04 
% against other code using 
% paths from outS.mat
%
CF = zeros(size(S));
CF(end,:) =max(K-S(end,:),0);
%
% Loop backwards in time 
%
for i = size(S)-1:-1:1
    CFidx = logical(max(K-S(i,:),0));

    X = S(i,CFidx)';               % Select in the money subset 
    Y = CF(i+1,CFidx)'*exp(-r*dt); % Y cash flow predictor
    
    %
    % Laguerre polynomial regression, the first 5 terms
    % Reference: the Longstaff-Schwartz paper
    % This part is also inspired from: Pricing American Options, by Mark Hoyle
    %
    Xn = X/S0; % Normalize paths
    A = [ exp(-Xn/2), (1-Xn).*exp(-Xn/2), (1-2*Xn+Xn.^2./2).*exp(-Xn/2), ...
          exp(-Xn/2).*1/6.*( -Xn.^3 + 9.*Xn.^2 - 18.*Xn + 6  )];
    %      exp(-Xn/2)./24.*( Xn.^4 - 16.*Xn.^3 + 72.*Xn.^2 - 96.*Xn + 24  ) ];
    b = A\Y;  
    CON = A*b; % Continuation 
    %
    %
    
%    pval  = polyfit(X,Y, 2); % Regress using second order polynomial
%    CON  = polyval(pval, X); % Continuation
 
    m_idx = max(K-X,0) > CON; % Test for continuation
        
    idx = false(1, N); % Form an empty logical set
    idx(CFidx) = m_idx;            % Populate it with indecies m_idx
    
    CF(i,idx) = max(K-X(m_idx),0);      % Set value to payoff otherwise
    CF(i,~idx) = exp(-r*dt)*CF(i+1,~idx);% push cashflow one step back
end
U(si) = mean(CF(1,:));%*exp(-r*dt);
end

U=U';
