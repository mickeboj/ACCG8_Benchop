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

function [xc,xtype]=GetNodesSpread(s2max,N2,s1max,N1,K)
% Points are generated in a triangle with corners (0,0) (0,s2max), (s2max,s1max)
% Domain goes from 0 to s2max and 0 to s1max
% Assuming N1 points at the top and N2 points vertically.
% Not working properly for K > 0 right now  
  fac = s1max/(K+s2max);
  x2 = linspace(0,s2max,N2);
  L0 = s1max;
  h0 = L0/(N1-1);
  xc = zeros(1,2);
  xtype(1)=1;
  for k=2:length(x2)
    n = ceil(fac*(K+x2(k))/h0 + 1);
    n = max(n,2);
    xx1 = linspace(0,fac*(K+x2(k)),n);
    xx2 = x2(k)*ones(1,n);
    xc = [xc; xx1(:) xx2(:)];
    xt = zeros(n,1);
    xt(1)=1;
    xt(end)=2;
    if (k==length(x2))
      xt(2:end-1) = 3;
    end
    xtype=[xtype;xt];
  end
  

