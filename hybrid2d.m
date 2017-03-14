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
function [disp]=hybrid2d(nodalcoord,nodalconn,fixeddofs,force,matprop,...
                         twodtype,nele,nnodes)
% This function performs the analysis for 2D linear finite 
% elements elasticity analysis with four node elements
% Input:
% an example mesh used to illustrate the organisation of the datatypes 
%
%    /|\  force 1 N 
%     |
%     |
%     |
%       
%  3-----4
%   | 2 |
%  2-----5
%   | 1 |
%  1-----6
% /\/\/\/\/\  fixed base
% 
% nodalcoord    - nodal coordinates in the form of  x y. The nodes are 
%                 assumed to be arranged in ascending order
% 
%     Example
%     1.000  2.000 % Node 1
%     1.000  4.000 % Node 2
%     1.000  6.000 % Node 3
%     2.000  6.000 % Node 4
%     2.000  4.000 % Node 5
%     2.000  2.000 % Node 6
% 
% nodalconn     - nodal connectivity in the anticlockwise fashion. The 
%                 elements are arranged to be in ascending order
% 
%     Example
%     1 6 5 2 % Element 1
%     2 5 4 3 % Element 2
% 
% 
% fixeddofs     - fixed degrees of freedom (with zero displacement)
%                 the degrees of freedom that are fixed
% 
%     Example 
%     [1 2 11 12] % Degrees of freedom associated with nodes 1,6
% 
% force         - force vector (make it sparse matrix)
% 
%     Example 
%     [0 0 0 0 0 0.5 0 0.5 0 0 0 0]'%dof associated with each node in order
% 
% matprop       - material properties in the form [E nu].
%                 Youngs Modulus and Poissons ratio
% 
%     Example     
%     [210e9 0.3]
%     
% twodtype      - the type of 2d analysis
%                 0  - plane strain
%                 1  - plane stress
% 
% nele          - total number of elements
% 
%     Example     
%     2
% 
% nnodes        - total number of nodes
% 
%     Example     
%     6
%     
% Output:
%
% disp          - The Displacement Vector



    %Initialize arrays

    %The Gauss Points
    c=1/sqrt(3);
    gauss = [c,  c; 
             c, -c;
            -c,  c;
            -c, -c];
   
    %Array for the Stress shape functions
    sshape=zeros(3,20);

    %Array for the Displacement shape functions
    dispshp=zeros(1,16);

    %Derivative of Displacement shape functions
    derdispshape=zeros(2,16);

    %Elemental degrees of freedom
    eledof=zeros(8,1);

    %Total number of degrees of freedom
    ndofs=2*nnodes;

    %Arrays to store the sparse array for the stiffness matrix
    %there will be duplication of values in the rows and columns
    %but matlab adds it up
    X=zeros(nele*64,1);Y=zeros(nele*64,1);Z=zeros(nele*64,1);

    %pointer to the next row to be filled up
    ntriplets=0;

    %The Displacement array
    disp=zeros(ndofs,1);

    %Form the Constitutive Matrices
        %Choice of 1 as the third argument selects plane stress
        %Choice of 0 as the third agrument selects plane strain
    [~,invC]=constitutive(matprop(1),matprop(2),twodtype);
    
    %Compute the values of the shape functions and their derivatives at the Gauss points
    for i = 1:4
            
            sshape(:,(i-1)*5+1:i*5)=stressshape(gauss(i,:));
            [dispshp(:,(i-1)*4+1:i*4),derdispshape(:,(i-1)*4+1:i*4)]=dispshape(gauss(i,:));
            
    end

    % Form and Assemble the Element Stiffness Matrices for every element

    for ele=1:nele
    
        %nodalconnectivity specific to this element
        nodconn=nodalconn(ele,:);
    
        %degrees of freedom beloning to this element
        eledof([1 3 5 7])=2*nodconn-1;    
        eledof([2 4 6 8])=2*nodconn;    
    
        % Call for the Element Stiffness matrix
        [Ke]=elestiff(invC,derdispshape,sshape,nodalcoord(nodconn,:));
    
        % Assembly
        for i = 1 : 8        
            for j = 1:8            
            
                ntriplets=ntriplets+1;            
                X(ntriplets)=eledof(i);            
                Y(ntriplets)=eledof(j);            
                Z(ntriplets)=Ke(i,j);            
                
            end        
        end
    
    end

    % Form the Sparse Matrix

    K=sparse(X,Y,Z,ndofs,ndofs);

    %Boundary Conditions
    alldofs=1:ndofs;
    freedofs=setdiff(alldofs,fixeddofs);

    %Solve for displacements
    disp(freedofs)=K(freedofs,freedofs)\force(freedofs);

