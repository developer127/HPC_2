clc, close all, clear all

% FEM2D  two-dimensional finite element method for linear second order PDE.
% Initialisation
load coordinates.dat; 
load elements.dat;
eval('load neumann.dat;','neumann=[];');
eval('load material.dat;','material=ones(size(elements,1),1);');
dirichlet = load('dirichletPoisson.dat');

%*** Refinement
for k=1:5
%[coordinates,elements,material,dirichlet,neumann] ...
           %= refineR(coordinates,elements,material,dirichlet,neumann);
[coordinates,elements,material,dirichlet] ...
           = refineR(coordinates,elements,material,dirichlet);
end
%*** Assembly of stiffness matrix
A = sparse(size(coordinates,1),size(coordinates,1));
for j = 1:size(elements,1)
  A(elements(j,:),elements(j,:)) = A(elements(j,:),elements(j,:)) ...
       + stima(coordinates(elements(j,:),:),material(j));
end
%*** Right hand side
b = zeros(size(coordinates,1),1);
for j = 1:size(elements,1)
  vertices = coordinates(elements(j,:),:);
  xM = sum(vertices)/3;                     %center of Gravity
  b(elements(j,:)) = b(elements(j,:)) ...
   +det([1 1 1;vertices']) * f1(xM,material(j))/6 ...
   -1/2*(vertices([2,3,1],:)-vertices([3,1,2],:))*[0,-1;1,0]*f2(xM,material(j))';
end
%*** Neumann conditions
for j = 1 : size(neumann,1)
   n = coordinates(neumann(j,2),:) - coordinates(neumann(j,1),:);
   dist = norm(n);
   n = ([0,1;-1,0]*n')/dist;
   xM = sum(coordinates(neumann(j,:),:))/2;  % Middle of the line
   b(neumann(j,:))=b(neumann(j,:)) + 0.5 * dist * (g(xM)+f2(xM,material(j))*n);
end
%*** Prescribe values at Dirichlet nodes
u = zeros(size(coordinates,1),1);
dirichlet = unique(dirichlet);
u(dirichlet) = uD(coordinates(dirichlet,:));
b = b - A * u;
%*** Computation of the solution
freenodes = setdiff(1:size(coordinates,1), dirichlet);
u(freenodes) = A(freenodes,freenodes) \ b(freenodes);
%*** Graphic representation
trisurf(elements,coordinates(:,1),coordinates(:,2),u,'facecolor','interp')

view(-16,20)
print -depsc2 'Solution'
title('Solution to L-Shape with mixed BC')

d_max = eigs(A(freenodes,freenodes))
d_min = eigs(A(freenodes,freenodes),5,'sm')

cond = d_max(5)/d_min(5)