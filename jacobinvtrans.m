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
function [J,DetJ,InvJ,Trans]=jacobinvtrans(derdispshape,nodalcoord)
%To form the Jacobian, its determinant, Inverse of the Jacobian
%and the Transformation Tensor for the stresses at a given point
%Input:
%   derdispshape   - The derivative of the displacement shape function at the
%                    point where the output terms are to be evaluated
%   nodalcoord     - The nodalcoordinates of the element in the order of the 
%                    nodal connectivity
%Output:
%   J              - Jacobian
%   DetJ           - Determinant of the Jacobian
%   InvJ           - The Inverse of the Jacobian
%   Trans          - The Transformation tensor for the stresses


    J=derdispshape*nodalcoord;    
    
    DetJ=det(J);

    InvJ=[J(2,2) -J(1,2);
          -J(2,1) J(1,1)]/DetJ;

    Trans=[(J').^2 2*(J(1,:).*J(2,:))';    
           (J(:,1).*J(:,2))' J(1,1)*J(2,2)+J(1,2)*J(2,1)];
       
end

