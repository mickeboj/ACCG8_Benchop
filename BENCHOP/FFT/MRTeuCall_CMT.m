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

function [U,errRel]=MRTeuCall_CMT(S,K,T,r,sig,lambda,gam,delta)
S=S(:)';
U=fftEU(@MertonKf,K,S,r,T,[sig lambda gam delta]);


function [P0]=fftEU(Kf,K0,S0,r,T,par)
% tic
m=2^8;
Zmax=2*pi/sqrt(par(1)*T)/40*m;  % replace sigma^2 in B&S with V
dz=Zmax/m;
z=(0:(m-1))*dz;
dx=2*pi/dz/m;
b=m/2*dx;
x=b-(0:(m-1))*dx;
S=exp(x+log(K0));
alpha=5;
w=ones(1,m);w([1 m])=1/2;
CC=exp(1i*b*z)./(1i*z+alpha)./(1i*z+alpha+1).*exp(Kf(z-1i*(alpha+1),T,r,0,par)-r*T);
P=dz*exp(alpha*x+x+log(K0)).*real(fft(CC.*w))/pi;
F=griddedInterpolant(fliplr(S),fliplr(P),'spline');
P0=F(S0);
% toc

% Compute relative error (compare with tabled results
% MRT_Exact=[7.542526012669795 14.336250885310660 22.359969565189871];
% errRel=(P0-MRT_Exact)./MRT_Exact

% Plot solution
% SGrid=60:1:160;
% plot(SGrid,F(SGrid));
% legend('Carr-Madan FFT')
% xlabel('S')
% ylabel('Price')

function [K,AI]=MertonKf(u,T,r,s0,par)
%function [K,AI]=MertonCf(u,T,r,s0,par)
%
% log characteristic function for the Merton model
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
% s=par(1);     : vol BM
% lam=par(2);   : jump int
% muj=par(3);   : jump mean
% sj=par(4);    : jump vol
%
% ---------------output-------------------
% 
%      K : log characteristic function 
%   anint: strip of analyticity for current parameter values 

% (C) 2006 Magnus Wiktorsson

s=par(1); % BM-vol
lam=par(2); % jump int
muj=par(3); % jump mean
sj=par(4);  % jump vol
mum=-lam*(exp(sj.^2/2+muj)-1);
K=(1i*u).*(s0+T.*(r+mum-s^2/2))+T.*(1i*u).^2/2*s^2+lam*T.*(exp((1i*u).^2/2*sj.^2+(1i*u)*muj)-1);
AI=[-1,1]/sj*10; %Due to numerical infinities
