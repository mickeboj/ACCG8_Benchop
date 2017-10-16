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

function [U,errRel]=HSTeuCall_CMT(S,K,T,r,V,kap,th,sig,rho)
S=S(:)';
U=fftEU(@HestonKf,K,S,r,T,[V  kap th sig rho]);

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
% HST_Exact=[2.302535842814927 7.379832496149447 14.974005277144057];
% errRel=(P0-HST_Exact)./HST_Exact


% % Plot solution
% SGrid=60:1:160;
% plot(SGrid,F(SGrid));
% legend('Carr-Madan FFT')
% xlabel('S')
% ylabel('Price')


function [K,anint,anint2]=HestonKf(z,T,r,s0,param)
%function [K,anint]=HestonKf(z,T,r,s0,param)
%
% log characteristic function for the Heston model
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
%  v0=param(1);    : Initial vol.
%  kap=param(2);   : Mean reversion
%  th=param(3);    : Mean level
%  s=param(4);     : Vol diffusion.
%  rho=param(5);   : Correlation
%
% ---------------output-------------------
% 
%      K : log characteristic function 
%   anint: strip of analyticity for current parameter values 

% (C) 2006 Magnus Wiktorsson
 
v0=param(1);    % Initial vol.
kap=param(2);   % Mean reversion
th=param(3);    % Mean level
s=param(4);     % Vol diffusion.
rho=param(5);   % Correlation

z2=-1i;
anint2=sort(roots([1 -(s^2-2*rho*s*kap)/(1-rho^2)/s^2 -kap^2/(1-rho^2)/s^2]))';
if nargout>=2
 y=max(T(:));
 sigma=s;
 kappa=kap;
 a2=roots([(-rho^3*sigma^3*y^3+rho*sigma^3*y^3) (3*kappa*rho^2*sigma^2*y^3-rho*sigma^3*y^3-kappa*sigma^2*y^3+6*rho^2*sigma^2*y^2-6*sigma^2*y^2) (-3*kappa^2*rho*sigma*y^3+kappa*sigma^2*y^3-12*kappa*rho*sigma*y^2+6*sigma^2*y^2-24*rho*sigma*y) kappa^3*y^3+6*kappa^2*y^2+24*kappa*y+48]);
 anint=[max(a2(a2<=0)) min(a2(a2>0))];
end 
d=sqrt((kap-rho*s*1i*z).^2+s^2*(1i*z-(1i*z).^2));
CC=kap*th/s^2*((kap-rho*s*1i*z).*T-2*log(((kap-1i*rho*s*z)+d.*coth(d.*T/2))./(d./sinh(d.*T/2))));
DD=((1i*z).^2-1i*z)./((kap-1i*rho*s*z)+d.*coth(d.*T/2));
K=1i*z.*(s0+r.*T)+CC+DD*v0;