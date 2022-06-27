function [eps_cr, eps_ii, sigma, x] = OneStepNortonArmstrongFrederick(eps, neps_cr, neps_ii, k, mu, A, n, kappa, c, dt)
sigma=k*trace(eps-neps_cr)*eye(3)+2*mu*dev(eps-neps_cr);%пробное напряжение
x=c*dev(neps_cr-neps_ii); %пробное микронапряжение
sigma_eff=sigma-x;
norm_stress_eff_dev=nor(dev(sigma_eff));
eps_cr=neps_cr+dt*A*norm_stress_eff_dev^(n-1)*dev(sigma_eff);
norm_eps_cr_dot=nor((eps_cr-neps_cr)/dt);%скорость деформации ползучести скаляр
eps_ii=neps_ii+dt*kappa*norm_eps_cr_dot*x;
sigma=k*trace(eps-eps_cr)*eye(3)+2*mu*dev(eps-eps_cr);%окончательное вычисление напряжений
end