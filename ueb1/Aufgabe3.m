clear all
% 1D poisson problem:

% Anzahl der Gitterpunkte n
k = [10:10:90, 100:50:250, 300:100:500];
cond_M = zeros(size(k));
cond_N = cond_M;

for i = 1:size(k,2)
    n = k(i);
    % M System Matrix f?r Dirichlet Randbedingungen
    M = spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);

    % N System Matrix f?r Direchlet am linken Rand und Neumann Randbedingungen
    %   am rechten Rand
    N = M;
    N(n,n-1)=-2;

    % Berechnung der Eingenwerte von M und N:
    % 1) Analytisch, 2) Numerisch:
    m_ana = 4*sin(pi/2*[1:n]'/(n+1)).^2;
    m_num = unique([eigs(M,5),eigs(M,5,'sm')]);
    cond_M(1,i) = m_ana(n)/m_ana(1);
    n_num = unique([eigs(N,5),eigs(N,5,'sm')]);
    cond_N(1,i) = n_num(10)/n_num(1);
end
hold on
title('Condition der Steifigkeitsmatritzen')
xlabel('Gitterpunkte')
ylabel('Condition Number')
plot(k,cond_M,'r','LineWidth',2)
plot(k,cond_N,'LineWidth',2)
legend('Dirichlet BC', 'Dirichlet BC left and Neumann BC right')
hold off