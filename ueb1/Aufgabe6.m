clear all
% 2D poisson problem:

% Anzahl der Gitterpunkte n
k = [10:10:100];
cond_A = zeros(size(k));

for i = 1:size(k,2)
    n = k(i);
    % T SubMatrix for Kronecker Product
    T = spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);
    
    % Assemblierung der Steifigkeitsmatrix A:
    A = kron(speye(n),T) + kron(T,speye(n));

    % Berechnung der Eingenwerte von M und N:
    % 1) Analytisch, 2) Numerisch:
    A_ana = 4*(sin(pi/2*[1:n]'/(n+1)).^2*ones(1,n)...
            + ones(n,1)*sin(pi/2*[1:n]/(n+1)).^2);
    A_ana = unique(A_ana(:));
    A_num = unique([eigs(A,5),eigs(A,5,'sm')]);
    cond_A(1,i) = A_ana(end)/A_ana(1);
end
hold on
title('Condition der Steifigkeitsmatrix')
xlabel('Gitterpunkte')
ylabel('Condition Number')
plot(k,cond_A,'LineWidth',2)
legend('Dirichlet BC')
hold off