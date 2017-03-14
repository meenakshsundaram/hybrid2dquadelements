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
%% A demo of a wrapper for a 2D Four Node Plane Stress or Strain Analysis of a Gripper using the hybrid formulation

%The final result
%ans =
%   99.8978
%  -99.8987
%
%
%ans =
%
%   81.8304
%   81.8317


%% Initializaiton
clc;
clear;
% Loading
% fixednodes    - nodes which need to be fixed
load demogripper/fixednodes.txt;   
% output nodes  - output nodes at which the displacements need to be found
load demogripper/outputnodes.txt;
% forcenodes    - nodes at which force needs to be applied
load demogripper/forcenodes.txt;
% forceval      - value of forces at the force nodes
load demogripper/forceval.txt
% matprop       - material properties
load demogripper/matprop.txt;
% nodalconn     - nodal connectivity
load demogripper/nodalconn.txt
% nodcalcoord   - nodal coordinate
load demogripper/nodalcoord.txt

% For further information of the pattern in which the data is organised for
% 1,2,3 look at the function file hybrid 2d

%% Analysis    

%the total number of elements
nele=size(nodalconn,1);

%the total number of nodes
nnodes=size(nodalcoord,1);

%the fixed degress of freedom from the fixed nodes
%fixeddofs contains the number of the degree of freedom which is fixed
fixeddofs=[2*fixednodes 2*fixednodes-1];

%the total number of degress of freedom
ndofs=2*nnodes;

%forming the force vector
force=sparse(ndofs,1);
force(2*forcenodes-1)=forceval;

%calling the analysis routines for hybrid formulation
[disp]=hybrid2d(nodalcoord,nodalconn,fixeddofs,force,matprop,1,nele,nnodes);

%post processing
h=figure(1);set(gca,'fontsize',32,'linewidth',2);
%plots the mesh and not very robust
plotmesh(nele,nodalcoord,1,disp,nodalconn,h,[1,0,0],[0,0,1]);

%displays the displacement of the output degrees of freedom
display(disp(2*outputnodes,:));
display(disp(2*outputnodes-1,:));
