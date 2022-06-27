function [res] = resid(ksi,eps, neps_cr, neps_ii, k, mu, A, n, kappa, c, dt)
F=(ksi/(A*dt))^(1/n);
F=max(F,1e-10);
eps_cr=UpdateEpscr(eps, neps_cr, neps_ii, ksi,kappa, c, mu, F);
eps_ii=UpdateEpsii(neps_ii, eps_cr,ksi,kappa,c);
sigma=k*trace(eps-eps_cr)*eye(3)+2*mu*dev(eps-eps_cr);%актуальное напряжение
x=c*dev(eps_cr-eps_ii); %актуальное микронапряжение
sigma_eff=sigma-x;
%наивный вариант
%res=(ksi/(A*dt))^(1/n)-nor(dev(sigma_eff));%невязка
%более эффективный вариант
%res=(ksi/(A*dt))-(nor(dev(sigma_eff))/(1-nomega))^n;%невязка
res=(ksi/(A*dt))-(nor(dev(sigma_eff)))^n;
end