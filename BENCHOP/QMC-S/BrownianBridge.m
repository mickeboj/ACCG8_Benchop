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
% Implementation of the Brownian bridge with midpoint division
% Code is based on the algorithmon p.85 in P. Glasserman, 
% Monte Carlo methods in Financial Engineering
% 
% By: Josef Hook, joh@kth.se
% Updated 2015, josef.hook@it.uu.se
%
% Usage: 
%
% Brownian bridge from a=-1 at time 0.1 to b=3 at time 5 in 8 steps
%
% figure,plot(BrownianBridge(-1,0.1,3,5, randn(1,2^3-1)))
%
% Ensemble example:
% figure,
% plot(BrownianBridge(zeros(40,1),0.1,sqrt(5).*randn(40,1),5, randn(40,2^8-1))')
% 
% Returns Wiener paths including Wt1 and Wt2
%
%
function [ret,t] = BrownianBridge(Wt1, t1, Wt2, t2, Z)


N = size(Z,1);
I = size(Z,2)+1; % Should always be odd since we need even intervalls


if(abs(rem(log2(I),1))>0)
    error('length of arg 5 should be on the form 2^m -1') 
end

%T = t2-t1;%t(end);
h = I;
jmax = 1;
t = linspace(t1, t2, I+1);
% output
w=zeros(N,I+1); % include left and right border
w(:,end) = Wt2;
w(:,1) = Wt1; 

ll=1;
for k = 1:log2(I)
    
    imin = h/2;
    i = imin;
    l = 0;
    r=h; 
    for j=1:jmax
        a = (   (t(r+1) - t(i+1)).*w(:,l+1) + (t(i+1) - t(l+1)).*w(:,r+1) )./(t(r+1) - t(l+1));
        b = sqrt( (t(i+1) - t(l+1)).*( t(r+1) - t(i+1) )./( t(r+1) - t(l+1) ) );
        
        w(:,i+1) = a + b.*Z(:,ll);
        
        ll=ll+1;
        i = i + h;
        l = l + h;
        r = r + h;
    
    end
    jmax = 2*jmax;
    h = imin;
    
end

ret = w; 