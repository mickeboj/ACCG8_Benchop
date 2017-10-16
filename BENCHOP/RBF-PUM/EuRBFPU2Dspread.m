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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.

function U=EuRBFPU2Dspread(S,K,T,r,sig1,sig2,rho,payoff,phi,ep,xc,M,npx,npy,op)

% Compute the value of a European spread option using uniform RBF-PUM approximation
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
% sig1 volatility on 1st asset x direction
% sig2 volatility on 2nd asset y direction
% rho correlation
% payoff payoff function
% op number of derivative
% npx number of partitions in x direction
% npy number of partitions in y direction
%
%  (C) Victor Shcherbakov & Elisabeth Larsson

% Define domain parameters
xeval = S; % evaluation points
Ne = size(S,1);

N = length(xc);

% Define boundaries
x_min = min(xc(:,1));
x_max = max(xc(:,1));
y_min = min(xc(:,1));
y_max = max(xc(:,2));
bxn = find(xc(:,1)==x_min); % near-field boundary in x
byn = find(xc(:,2)==y_min); % near-field boundary in y

h = (xc(3,2)-xc(2,2)); 

% Define weights on BDF2
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);

% Compute the size of the partitions
Hx = (x_max-x_min)/npx;
x = Hx/2:Hx:x_max-Hx/2;
Hy = (y_max-y_min)/npy;
y = Hy/2:Hy:y_max-Hy/2;
% overlap = 1.5*h; % overlap size
overlap = 5*h;

dimenx = min(diff(x));
dimeny = min(diff(y));
dimen = max(dimenx,dimeny);
R = hypot(0.5*dimen,0.5*dimen) + overlap; % radius of partitions

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
    
    B = (r*X*Dx + r*Y*Dy ...
        + 0.5*sig1^2*X.^2*Dxx ...
        + rho*sig1*sig2*X*Y*Dxy ...
        + 0.5*sig2^2*Y.^2*Dyy - r*D0)/A0;
  
    ind_e = iboxeval(i).indlocpts;
    flagempt = isempty(ind_e);
    if flagempt == 0
        re = xcdist(xeval(ind_e,:),xloc,op);
        Wei = spdiags(We{i},0,Nei(i),Nei(i));
        E{i} = Wei*RBFmat(phi,ep,re,op)/A0; % evaluation matrix
    end
    
    L(ibox(i).indlocpts,ibox(i).indlocpts) = L(ibox(i).indlocpts,ibox(i).indlocpts) + B;
   
end

% Zero in the rows corresponding to boundary conditions
L(bxn,:) = 0;
L(byn,:) = 0;

% Time-stepping matrix I-beta_0*L in the interior
C = speye(N) - beta0*L;

[Lc,Uc] = lu(C);

tn = k(1);
u0 = max(xc(:,1)-xc(:,2),0);
size(u0);

rhs = u0;
rhs(bxn) = u0(bxn);
rhs(byn) = u0(byn);

% Time-stepping loop
for n = 1:M
    
    % New solution value from linear system
    u = Uc\(Lc\rhs);

    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    
    % boundary conditions
    rhs(bxn) = 0;
    rhs(byn) = xc(byn,1);
    
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










