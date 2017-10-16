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

function U=HestonRBF2D(S,K,T,r,V,kappa,theta,sig,rho,payoff,phi,ep,xc,xtype,M,op)
% Compute the value of a European option using uniform RBF approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
  
% The parameters are given unscaled, but in the solver we always scale
% down x and U so that K=1.
  N = length(xc);
  b0 = find(xtype==1);
  bx = find(xtype==2);
  Nv = length(b0);
  Ns = N/Nv;
  %
  % Define the scaling to get a certain ratio between the grid sizes in s
  % and v
  %
  K0 = K;
  hs = max(xc(:,1))/K0/(Ns-1);
  hv = 1/2.5*hs; % From tests, may be problem dependent
  % We want max(xc(:,2))/V0/(Nv-1) = hv
  V0 = max(xc(:,2))/(Nv-1)/hv;
  K = 1;
  %
  % Evaluation grid points
  %
  [ss,vv]=meshgrid(S/K0,V/V0);
  xe = [ss(:),vv(:)];
  xc(:,1) = xc(:,1)/K0;
  xc(:,2) = xc(:,2)/V0;
  %
  % The boundary is now 2-D and needs to be identified
  %
  x_max = max(xc(:,1));

  fac = K0; % Due to scaled payoff.


  [k,beta0,beta1,beta2]=BDF2coeffs(T,M);
  %
  % Diagonal matrices with coordinates for coefficients
  %
  X = spdiags(xc(:,1),0,N,N);
  Y = spdiags(xc(:,2),0,N,N);
  %
  % Distance matrix to build interpolation and collocation matrices
  %
  rc = xcdist(xc,xc,1);
  re = xcdist(xe,xc,0);
  %
  % Interpolation and operator matrices
  %
  A0 = RBFmat(phi,ep,rc,'0');
  Ax = RBFmat(phi,ep,rc,'1',1);
  Axx = RBFmat(phi,ep,rc,'2',1);
  Ay = RBFmat(phi,ep,rc,'1',2);
  Ayy = RBFmat(phi,ep,rc,'2',2);
  Axy = RBFmat(phi,ep,rc,'m2',[1 2]);
  Ae = RBFmat(phi,ep,re,op);
  %
  % Build the Black-Scholes operator matrix
  %
  I = speye(N);
  L = (0.5*V0*X.^2.*Y*Axx + rho*sig*X.*Y*Axy + sig^2/2/V0*Y*Ayy + r*X*Ax ...
       + kappa*(theta/V0*I-Y)*Ay - r*A0)/A0;
  
  clear Ax Axx Axy Ayy Ay
  %
  % The evaluation matrix A*A^-1
  % 
  Ae = Ae/A0;
  %
  % Time-stepping matrix I-beta_0*L in the interior
  %
  C = eye(N) - beta0*L;
  %
  % Unit one in the rows corresponding to Dirichlet boundary conditions
  %
  C(bx,:) = I(bx,:);
  C(b0,:) = I(b0,:);

  [Lc,Uc] = lu(C);
  %
  % Initial condition. Note that boundaries are at the next time
  %
  tn = k(1);
  u0 = payoff(xc(:,1),K,r,0);

  rhs(:,1) = u0;
  rhs(bx,1) = payoff(x_max,K,r,tn);
  rhs(b0,1) = 0;
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
    rhs(bx) = payoff(x_max,K,r,tn);
    rhs(b0) = 0;
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

