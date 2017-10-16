 %Copyright (C) 2015 Erik Lindström

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

function [U,errRel]=BSeuCallDeltaII_CMT(S,K,T,r,sigma)
S=S(:)';
U=fftEU(@BSKf,K,S,r,T,sigma);

function [P0]=fftEU(Kf,K0,S0,r,T,par)
% tic, 
m=2^8;
Zmax=2*pi/sqrt(par^2*T)/40*m;
dz=Zmax/m;
z=(0:(m-1))*dz;
dx=2*pi/dz/m;
b=m/2*dx;
x=b-(0:(m-1))*dx;
S=exp(x+log(K0));
alpha=150;
w=ones(1,m);w([1 m])=1/2;
CC=exp(1i*b*z)./(1i*z+alpha).*exp(Kf(z-1i*(alpha+1),T,r,0,par)-r*T);
P=dz*exp(alpha*x).*real(fft(CC.*w))/pi;
F=griddedInterpolant(fliplr(S),fliplr(P),'spline');
P0=F(S0);
% toc

% % Compute relative error
% [~,D0]=BSCall_exact(S0,K0,T,r,par);
% errRel=(P0-D0)./D0

% % Plot solution
% SGrid=60:1:160;
% [~,DeltaGrid]=BSCall_exact(SGrid,K0,T,r,par);
% plot(SGrid,F(SGrid),SGrid,DeltaGrid,'--');
% legend('Carr-Madan FFT','Exact')
% xlabel('S')
% ylabel('Delta')

function [K,anint]=BSKf(u,T,r,s0,par)
% log characteristic function for the Black-Scholes model
%
% ------------input-----------------------
%
%     z: evaluation point, Note that we allow z to be complex, if we
%        evaluate at the point -ix (for x real) we evaluate the log
%        moment generating function at the point x. 
%     T: time
%     r: continously compounded short rate 
%    s0: intial/current log stock price 
% param: model parameters
%
% ---------------output-------------------
% 
%      K : log characteristic function 
%   anint: strip of analyticity for current parameter values 

% (C) 2006 Magnus Wiktorsson

s=par(1); % BM-vol
K=1i*u.*(s0+T.*(r-s.^2/2))-T.*u.^2/2*s.^2;
anint=[-inf,inf];
