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
function [sshape]=stressshape(coord)
%This function evaluates the stress shape function at a given coordinate
%Input:
%       coord  - Coordinate at which the stress shape function 
%                 is evaluated
%Output:
%       sshape - Stress Shape function

    sshape=zeros(3,5);    
    
    sshape(1,[1,2])=[1 coord(2)];    
    sshape(2,[3,4])=[1 coord(1)];    
    sshape(3,5)=1;