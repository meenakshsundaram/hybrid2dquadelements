%
%    This file is a part of the hybrid2dquadelements a matlab library. This 
%    is free software: you can redistribute it and/or modify it under 
%    the terms of the GNU Lesser General Public License as published by the 
%    Free Software Foundation, either version 3 of the License, or 
%    (at your option) any later version.
%
%    hybrid2dquadelements is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with hybrid2dquadelements.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2010 Meenakshi Sundaram   
%
function [C,invC]=constitutive(E,nu,flag)
% This function forms the constitutive tensor and its inverse
% 
% Input :
%     E     - Young's Modulus 
%     nu    - Poissons Ration
%     flag  - 1 indicates Plane Stress 
%             0 indicates Plane Strain
% Output:
%     C     - Constitutive Tensor
%     invC  - Inverse of the Constitutive Tensor
%

    if(flag)   
    
        %Plane Stress
        C=E/(1.0-nu^2)*[ 1.0    nu   0.0;
                         nu     1.0  0.0;
                         0.0    0.0  (1.0-nu)/2.0 ];            
                  
        invC=(1.0/E)*[ 1.0  -nu   0.0;
                      -nu    1.0  0.0;
                       0.0   0.0  2.0*(1.0+nu) ];
             
    else
    
        %Plane Strain
        invC=((1.0+nu)/E)*[ 1.0-nu  -nu      0.0;
                           -nu       1.0-nu  0.0;
                            0.0      0.0     2.0 ];
                 
         C= (E/((1.0+nu)*(1.0-2.0*nu)))*[ 1.0-nu  nu      0.0;
                                          nu      1.0-nu  0.0;
                                          0.0     0.0     (1.0-2.0*nu)/2.0 ];
                            
    end