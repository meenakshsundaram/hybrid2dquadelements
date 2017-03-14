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
function [dispshape,derdispshape]=dispshape(coord)
%This function evaluates the displacement shape functions and its
%derivative at a given coordinate
%Input:
%   coord          - The coordinate at which the displacement shape function is
%                    evaluated
%Output:
%   dispshape      - The displacement shape function 
%   derdispshape   - The derivative of the displacement shape function


    RP = 1.0+coord(1);
    SP = 1.0+coord(2);
    RM = 1.0-coord(1);
    SM = 1.0-coord(2);
    
    dispshape = 0.25*[ RP*SP RM*SP RM*SM RP*SM ];
    
    derdispshape = 0.25*[ SP -SP  -SM  SM;
                          RP  RM  -RM -RP ];
                     
