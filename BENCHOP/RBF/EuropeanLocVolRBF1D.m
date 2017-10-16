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

function U=EuropeanLocVolRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep,limited);
  % Compute the value of a European option using uniform RBF approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
  
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
  %
  % Create a modified x-vector to use with the sigma function if needed
  % This is for the case where the sigma function is not valid for small s
  %
  xsig = xc;
  if (strcmpi(limited,'yes'))
    xsig = max(xsig,0.05*K);
  end
 
  x_min = 0;

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
  % The differentation matrix A*A^-1
  % 
  Ae = Ae/A0;
  %
  % Initial condition. Note that boundaries are at the next time
  % 
  u0 = payoff(xc,K,r,0);  
  rhs(:,1) = u0;
  tn = k(1);
  rhs(end,1) = payoff(x_max,K,r,tn);
  %
  % Build the operator matrices that we will need
  %
  Lxx = 0.5*X.^2*Axx/A0;
  Lx0 = r*(X*Ax/A0  - eye(N));
  %
  % Time-stepping loop
  %
  C = zeros(N,N);
  C(1,:) = 0;
  C(end,:) = 0;
  C(1,1) = 1;
  C(end,end) = 1;
  I = speye(N);
  for n = 1:M
    %
    % Build the Black-Scholes operator matrix for the current time step
    %
    Sigma = spdiags(sig(xsig*K0,max(0,T-tn)),0,N,N);
    L = Sigma.^2*Lxx + Lx0;
    %
    %
    % Time-stepping matrix I-beta_0*L in the interior
    %
    C(2:N-1,:) = I(2:N-1,:) - beta0*L(2:N-1,:);
    %
    % New solution value from linear system
    %
    u = C\rhs;
    %
    % Prepare the right hand side for the next step
    %
    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    tn = max(0,tn);
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    % boundary conditions
    rhs(1) = 0;
    rhs(end) = payoff(x_max,K,r,tn);
    u0 = u;
  end
  %
  % Compute the solution or some derivative at the evaluation points
  %
  U = Ae*u;
  %
  % Compare with the reference values
  %
  U = U*fac;

