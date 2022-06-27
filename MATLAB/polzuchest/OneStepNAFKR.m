function [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKR(eps, neps_cr, neps_ii, nomega, k0, mu0, A0, n_Norton, kappa0, c0, B, m, dt)
%явная схема для модели с повреждаемостью
sigma=(1-nomega)*(k0*trace(eps-neps_cr)*eye(3)+2*mu0*dev(eps-neps_cr));%пробное напряжение
x=(1-nomega)*c0*dev(neps_cr-neps_ii); %пробное микронапряжение
sigma_eff=sigma-x;
norm_stress_eff_dev=nor(dev(sigma_eff));
eps_cr=neps_cr+(1-nomega)^(-n_Norton)*dt*A0*norm_stress_eff_dev^(n_Norton-1)*dev(sigma_eff);
norm_eps_cr_dot=nor((eps_cr-neps_cr)/dt);%скорость деформации ползучести скаляр
eps_ii=neps_ii+(1-nomega)^(-1)*dt*kappa0*norm_eps_cr_dot*x;
omega=nomega+dt*B*(1-nomega)^(-m)*norm_stress_eff_dev^m;
sigma=(1-omega)*(k0*trace(eps-eps_cr)*eye(3)+2*mu0*dev(eps-eps_cr));%окончательное вычисление напряжений
x=(1-omega)*c0*dev(eps_cr-eps_ii); %окончательное микронапряжение
end