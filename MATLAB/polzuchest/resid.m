function [res] = resid(ksi,eps, neps_cr, neps_ii, k, mu, A, n, kappa, c, dt)
F=(ksi/(A*dt))^(1/n);
F=max(F,1e-10);
eps_cr=UpdateEpscr(eps, neps_cr, neps_ii, ksi,kappa, c, mu, F);
eps_ii=UpdateEpsii(neps_ii, eps_cr,ksi,kappa,c);
sigma=k*trace(eps-eps_cr)*eye(3)+2*mu*dev(eps-eps_cr);%���������� ����������
x=c*dev(eps_cr-eps_ii); %���������� ���������������
sigma_eff=sigma-x;
%������� �������
%res=(ksi/(A*dt))^(1/n)-nor(dev(sigma_eff));%�������
%����� ����������� �������
%res=(ksi/(A*dt))-(nor(dev(sigma_eff))/(1-nomega))^n;%�������
res=(ksi/(A*dt))-(nor(dev(sigma_eff)))^n;
end