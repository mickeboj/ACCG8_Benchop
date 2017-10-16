 %Copyright (C) 2015 Elisabeth Larsson

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

function Phi=EuropeanVega(S,K,T,r,sig,payoff,phi,ep,xc,M,varep)
% We will compute u and vega simultaneously, since the right hand side
% for Vega depends on the second derivative of u.
% Compute the value of a European option using uniform RBF approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
  
  op = '0';
% The parameters are given unscaled, but in the solver we always scale
% down x and U so that K=1.
  K0 = K;
  K = 1;
  %
  % Evaluation grid points (uniform)
  %
  xe = S/K0;
  xc = xc/K0;
  N = length(xc);
  x_max = xc(end);
  p = str2num(op);
  fac = K0^(1-p);

  x_min = 0;
  sigma2 = sig^2;

  [k,beta0,beta1,beta2]=BDF2coeffs(T,M);
  %
  % If variable epsilon is requested
  %
  if strcmp(lower(varep),'yes')
    d(1)=xc(2)-xc(1);
    d(2:N-1,1) = 0.5*(xc(3:N)-xc(1:N-2));
    d(N) = xc(N)-xc(N-1);
    d = d/max(d);
    % We let the scaling be between ep/Q and ep
    % ep(dmax)=ep/Q, ep/Q/(dmin+alpha)=ep, dmin+alpha = 1/Q
    Q=16;
    ep0=ep;
    alpha = max(0,1/Q-min(d));
    ep = ep0/Q./(alpha+d);
  end  
  %
  % Diagonal matrices with coordinates for coefficients
  %
  X = spdiags(xc(:,1),0,N,N);
  %
  % Distance matrix to build interpolation and collocation matrices
  %
  rc = xcdist(xc,xc,1);
  re = xcdist(xe,xc,p);
  %
  % Interpolation and operator matrices
  %
  A0 = RBFmat(phi,ep,rc,'0');
  Ax = RBFmat(phi,ep,rc,'1',1);
  Axx = RBFmat(phi,ep,rc,'2',1);
  Ae = RBFmat(phi,ep,re,op);
  %
  % Build the Black-Scholes operator matrix
  %
  L = (sigma2/2*X.^2*Axx + r*X*Ax  - r*A0)/A0;
  %
  clear Ax 
  %
  % The differentation matrix A*A^-1
  % 
  Ae = Ae/A0;
  Axx = Axx/A0;
  %
  % Time-stepping matrix I-beta_0*L in the interior
  %
  C = eye(N) - beta0*L;
  %
  % Unit one in the rows corresponding to boundary conditions
  %
  C(1,:) = 0;
  C(end,:) = 0;
  C(1,1) = 1;
  C(end,end) = 1;
  [Lc,Uc] = lu(C);
  %
  % Initial condition. Note that boundaries are at the next time
  % 
  u0 = payoff(xc,K,r,0); 
  %
  % Need to think of asymptotic boundary conditions for phi, but probably zero.
  % Can check analytical formula.
  %
  phi0 = zeros(size(u0));
  phirhs = phi0;
  rhs(:,1) = u0;
  tn = k(1);
  rhs(end,1) = payoff(x_max,K,r,tn);
  %
  % Time-stepping loop
  %
  for n = 1:M
    %
    % New solution value from linear system
    %
    u = Uc\(Lc\rhs);
    %
    % Prepare the right hand side for the next step
    %
    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    % boundary conditions
    %    dudt = K0*(u-u0)/k(nextstep-1);
    %bsop = L*u;
    %force = K0*sig*X.^2*Axx*u;
    
    %figure(3),hold on,plot(xc*K0,dudt,xc*K0,K0*bsop,xc*K0,force)
    
    rhs(1) = 0;
    rhs(end) = payoff(x_max,K,r,tn);
    u0 = u;
    %
    % The computation of phi is one step after. We use the newly computed u
    % Think about the times. Perhaps write it down if it seems to work.
    % 
    phirhs = phirhs + beta0*sig*K0*X.^2*(Axx*u);
    phirhs(1)=0;
    phirhs(end)=0;
    phi = Uc\(Lc\phirhs);
    phirhs = beta1(nextstep)*phi-beta2(nextstep)*phi0;
    phi0 = phi;
  end
  %
  % Compute the solution or some derivative at the evaluation points
  %
  Phi = Ae*phi;
  %
  % What is the factor in this case? The same? Check. The right hand side
  % will be off by a factor.
  %
  %U = Ae*u*K0;
  %
  % Compare with the reference values
  %
  %U = U*fac;

