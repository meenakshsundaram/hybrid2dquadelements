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
function plotmesh(nele,nodedata,scale,disp,eledata,h,orig_color,disp_color)
%function plots the deformed and undeformed meshes
%Input:
%   nele        - number of elements
%   nodedata    - nodalcoordinates
%   scale       - the scaling for displacements
%   disp        - the displacements
%   eledata     - the element connectivity data
%   h           - the handle for the figure where 
%                 the plotting needs to be done
%   orig_color  - the color for the undeformed mesh
%   disp_color  - the color for the deformed mesh


    set(0,'CurrentFigure',h); hold on;
    for ele=1:nele
        
        elenodes=eledata(ele,1:4);
        h1=plot(nodedata(elenodes([1:4,1]),1),nodedata(elenodes([1:4,1]),2),'-','linewidth',3,'color',orig_color);        
        
    end
    
    nodedata=nodedata+[disp(1:2:end) disp(2:2:end)]*scale;   
    
    for ele=1:nele
        
        elenodes=eledata(ele,1:4);
        h2=plot(nodedata(elenodes([1:4,1]),1),nodedata(elenodes([1:4,1]),2),'-','linewidth',3,'color',disp_color);        
        
    end
    
    min_nodedata = min(nodedata);
    max_nodedata = max(nodedata);
    
    range_x = abs(max_nodedata(1)-min_nodedata(1));
    range_y = abs(max_nodedata(2)-min_nodedata(2));
        
    legend([h1,h2],'original mesh','deformed mesh','location','best');    
    axis([min_nodedata(1)-0.1*range_x max_nodedata(1)+0.1*range_x min_nodedata(2)-0.1*range_y max_nodedata(2)+0.1*range_y]);
    axis equal;
