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
function [Ke]=elestiff(invC,derdispshape,sshape,nodalcoord)
%This function forms the Element stiffness matrix
%Input :
% invC          - Inverse of the Constitutive Tensor
% derdispshape  - The Derivative of the Displacement Shape function at all
%                 Gauss Points
% sshape        - The Stress shape function at all Gauss Points
% nodalcoord    - The nodalcoordinates for the element in the order of
%                 connectivity
%Output:
% Ke            - Element Stiffness Matrix

    %Initialization

    B=zeros(3,8);
    G=zeros(5,8);
    H=zeros(5);

    % For every one of the gaussian points
    for i =1 :4
    
        %The derivative of the displacement is formed    
        derdisp=derdispshape(:,(i-1)*4+1:i*4);
    
        %Form the Jacobian, Determinant of Jacobian, Inverse of the Jacobian 
        %and the Transformation Tensor at each of the Gauss Points
        [~,DetJ,InvJ,Trans]=jacobinvtrans(derdisp,nodalcoord);
    
    
        if(DetJ<=0)
            display('error: negative jacobian');
        end
    
        %Formation of the Strain Displacement matrix
        B(1,[1 3 5 7]) = InvJ(1,:)*derdisp;
        B(2,[2 4 6 8]) = InvJ(2,:)*derdisp;
        B(3,[1 3 5 7]) = B(2,[2 4 6 8]);
        B(3,[2 4 6 8]) = B(1,[1 3 5 7]);
        
        %Formation of the P matrix
        P=Trans*sshape(:,(i-1)*5+1:i*5);
    
        %Accumulate G and H
        G=G+P'*B*DetJ;
        H=H+P'*invC*P*DetJ; 

    end

    %Find the Element Stiffness matrix
    Ke=G'*(H\G);

