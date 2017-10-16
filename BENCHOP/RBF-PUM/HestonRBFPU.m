 %Copyright (C) 2015 Victor Shcherbakov

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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.% What is the boundary condition for the spread option?
function U = HestonRBFPU(S,K,T,r,V,kappa,theta,sig,rho,payoff,phi,ep,xc,M,npx,npy,op,xtype)

% Compute the value of a European option under the Heston model
% using uniform RBF-PUM approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in x direction
% y_max is the size of the computational domain in y direction
% xc is computational points
% S evaluation points, points we want to know the price at
% K strike
% T maturity time
% r risk free interest rate
% V volatility 
% kap mean reversion
% th mean level
% sig vol diffusion
% rho correlation
% payoff payoff function
% op number of derivative
% npx number of partitions in x direction
% npy number of partitions in y direction
% xtype type of nodes (internal, boundary etc)
%  
% The parameters are given unscaled, but in the solver we always scale
% down x and U so that K=1.
%
%  (C) Victor Shcherbakov & Elisabeth Larsson

K0 = K;
K = 1;
  
[ss,vv]=meshgrid(S/K0,V);
xeval = [ss(:),vv(:)];
Ne = length(xeval);

xc(:,1) = xc(:,1)/K0;
xc(:,2) = xc(:,2);
N = length(xc);

h = max(diff(xc(:,1)));

% The boundary is now 2-D and needs to be identified
x_min = min(xc(:,1));
x_max = max(xc(:,1));
y_min = min(xc(:,2));
y_max = max(xc(:,2));
b0 = find(xtype==1);
bx = find(xtype==2);
by = find(xtype==3);
p = str2num(op);
fac = K0^(1-p);

% Coefficients for BDF-2
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);


% Compute the size of the partitions
Hx = (x_max-x_min)/npx;
x = Hx/2:Hx:x_max-Hx/2;
Hy = (y_max-y_min)/npy;
y = Hy/2:Hy:y_max-Hy/2;
overlap = 5*h; %1.5*h; % overlap size

dimenx = max(diff(x));
dimeny = max(diff(y));
dimen = max(dimenx,dimeny);
R = x_max/14;%hypot(0.5*dimen,0.5*dimen) + overlap; % radius of partitions x_max/8;% x_max/8;%

[xm,ym] = meshgrid(x,y); % centers of partitions
xm = [xm(:) ym(:)];
Np = length(xm); % number of partitions


% Check which point goes to which partition
for i=1:Np
    flagin = hypot(xm(i,1)-xc(:,1),xm(i,2)-xc(:,2)) < R;
    ibox(i).indlocpts = find(flagin);
    Ni(i)=length(ibox(i).indlocpts);
end



% Compute the weight functions for each partition using Shepard's method
s=zeros(N,1);
sx=zeros(N,1);
sy=zeros(N,1);
sxx=zeros(N,1);
sxy=zeros(N,1);
syy=zeros(N,1);

for i=1:Np
    [psi{i},psix{i},psiy{i},psixx{i},psixy{i},psiyx{i},psiyy{i}] = RBFPUweight2D(xc(ibox(i).indlocpts,:),xm(i,:),R);
    
    s(ibox(i).indlocpts) = s(ibox(i).indlocpts) + psi{i};
    sx(ibox(i).indlocpts)= sx(ibox(i).indlocpts) + psix{i};
    sy(ibox(i).indlocpts)= sy(ibox(i).indlocpts) + psiy{i};
    sxx(ibox(i).indlocpts) = sxx(ibox(i).indlocpts) + psixx{i};
    sxy(ibox(i).indlocpts) = sxy(ibox(i).indlocpts) + psixy{i};
    syy(ibox(i).indlocpts) = syy(ibox(i).indlocpts) + psiyy{i}; 
end  

% Define partition of unity
for i=1:Np
    s1 = s(ibox(i).indlocpts);
    s2 = s1.*s1;
    s3 = s2.*s1;
    sx1 = sx(ibox(i).indlocpts);
    sx2 = sx1.*sx1;
    sxx1 = sxx(ibox(i).indlocpts);
    sy1 = sy(ibox(i).indlocpts);
    sy2 = sy1.*sy1;
    syy1 = syy(ibox(i).indlocpts);
    sxy1 = sxy(ibox(i).indlocpts);


    W{i} =  psi{i}./s1;
    Wx{i} = psix{i}./s1 - psi{i}.*sx1./s2;
    Wy{i} = psiy{i}./s1 - psi{i}.*sy1./s2;
    Wxx{i} = -2*psix{i}.*sx1./s2 + psixx{i}./s1 ...
      + psi{i}.*(2*sx2./s3 - sxx1./s2);
    Wxy{i} = psixy{i}./s1 - psix{i}.*sy1./s2 ... 
      - psiy{i}.*sx1./s2 ...
      - psi{i}.*sxy1./s2 ... 
      + 2*psi{i}.*sx1.*sy1./s3;
    Wyy{i} = -2*psiy{i}.*sy1./s2 + psiyy{i}./s1 ... 
       + psi{i}.*(2*sy2./s3 - syy1./s2);
end


% Do the same for eval points
for i = 1:Np
    flagin = hypot(xm(i,1)-xeval(:,1),xm(i,2)-xeval(:,2)) < R;
    iboxeval(i).indlocpts = find(flagin);
    Nei(i)=length(iboxeval(i).indlocpts);
end

s=zeros(Ne,1);
for i = 1:Np
    flagempty = isempty(iboxeval(i).indlocpts);
    if flagempty == 0
        [psi{i},psix{i},psiy{i},psixx{i},psixy{i},psiyx{i},psiyy{i}] = RBFPUweight2D(xeval(iboxeval(i).indlocpts,:),xm,R);
        s(iboxeval(i).indlocpts) = s(iboxeval(i).indlocpts) + psi{i};
    end
end 

for i = 1:Np
    flagempty = isempty(iboxeval(i).indlocpts);
    if flagempty == 0
        We{i} =  psi{i}./s(iboxeval(i).indlocpts);
    end
end

% Form local RBF matrices
L = spalloc(N,N,sum(Ni.^2));
I = speye(N);
for i=1:Np
    
    xloc = xc(ibox(i).indlocpts,:);
    rc = xcdist(xloc,xloc,1);
    A0  = RBFmat(phi,ep,rc,'0');
    Ax  = RBFmat(phi,ep,rc,'1',1);
    Ay  = RBFmat(phi,ep,rc,'1',2);
    Axx = RBFmat(phi,ep,rc,'2',1);
    Axy = RBFmat(phi,ep,rc,'m2',[1,2]);
    Ayy = RBFmat(phi,ep,rc,'2',2);
    
    Wi = spdiags(W{i},0,Ni(i),Ni(i));
    Wix = spdiags(Wx{i},0,Ni(i),Ni(i));
    Wiy = spdiags(Wy{i},0,Ni(i),Ni(i));
    Wixx = spdiags(Wxx{i},0,Ni(i),Ni(i));
    Wixy = spdiags(Wxy{i},0,Ni(i),Ni(i));
    Wiyy = spdiags(Wyy{i},0,Ni(i),Ni(i));
    
    X = spdiags(xloc(:,1),0,Ni(i),Ni(i));
    Y = spdiags(xloc(:,2),0,Ni(i),Ni(i));
    
    D0 = Wi*A0; 
    Dx = Wi*Ax + Wix*A0;
    Dy = Wi*Ay + Wiy*A0;

    Dxx = Wi*Axx + 2*Wix*Ax + Wixx*A0; 
    Dyy = Wi*Ayy + 2*Wiy*Ay + Wiyy*A0;

    Dxy = Wi*Axy + Wix*Ay + Wiy*Ax + Wixy*A0;
    
    
    B = (0.5*X.^2.*Y*Dxx + rho*sig*X.*Y*Dxy + sig^2/2*Y*Dyy + r*X*Dx ...
       + kappa*(theta*I(ibox(i).indlocpts,ibox(i).indlocpts)-Y)*Dy - r*D0)/A0;
  
    ind_e = iboxeval(i).indlocpts;
    flagempt = isempty(ind_e);
    if flagempt == 0
        re = xcdist(xeval(ind_e,:),xloc,op);
        Wei = spdiags(We{i},0,Nei(i),Nei(i));
        E{i} = Wei*RBFmat(phi,ep,re,op)/A0; % evaluation matrix
    end
    
    L(ibox(i).indlocpts,ibox(i).indlocpts) = L(ibox(i).indlocpts,ibox(i).indlocpts) + B;
   
end
 
  
C = eye(N) - beta0*L;
C(bx,:) = I(bx,:);
C(b0,:) = I(b0,:);

[Lc,Uc] = lu(C);

tn = k(1);
u0 = payoff(xc(:,1),K,r,0);


rhs(:,1) = u0;
rhs(bx,1) = payoff(x_max,K,r,tn);
rhs(b0,1) = 0;

% Time-stepping loop
for n = 1:M
    
    % New solution value from linear system
    u = Uc\(Lc\rhs);

    % Prepare the right hand side for the next step
    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    
    % boundary conditions
    rhs(bx) = payoff(x_max,K,r,tn);
    rhs(b0) = 0;
    
    % update
    u0 = u;
end

% Compute the solution or some derivative at the evaluation points
U = zeros(Ne,1);
for i=1:Np
    ind=ibox(i).indlocpts;
    ind_e = iboxeval(i).indlocpts;
    flagempty = isempty(ind_e);
    if flagempty == 0
        U(ind_e) = U(ind_e) + E{i}*u(ind);
    end
end  

% Scale back
U = U*fac;

