clear all
% 1D poisson problem:

% Anzahl der Gitterpunkte n
n = 5;

% M System Matrix für Dirichlet Randbedingungen
M = spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);

% N System Matrix für Direchlet am linken Rand und Neumann Randbedingungen
%   am rechten Rand
N = M;
N(n,n-1)=-2;

% Berechnung der Eingenwerte von M und N:
% 1) Analytisch, 2) Numerisch:
m_ana = ones(n,1);
for k=1:n
    m_ana(k,1)= 4*sin(pi/2*k/(n+1))^2;
end
m_num = eigs(M);
m_dif = abs(m_num - m_ana)

n_ana = ones(n,1);
for k=1:n
    n_ana(k,1)= 4*sin(pi/2*k/n)^2;
end
n_num = eigs(N);
n_dif = abs(n_num - n_ana)