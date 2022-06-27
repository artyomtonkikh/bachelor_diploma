function [eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega] = OneStep(eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, mater_param, dt)
%input: eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, mater_param, dt
%output: eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega
d=sqrt(nomega^2+nf^2);
sigma=(1-d)*(mater_param.k*trace(eps-neps_cr-neps_p)*eye(3)+2*mater_param.mu*dev(eps-neps_cr-neps_p));
sigma_m=1/3*trace(sigma);
x1=(1-d)*mater_param.c1*dev(neps_cr+neps_p-neps_ii1);
x2=(1-d)*mater_param.c2*dev(neps_cr+neps_p-neps_ii2);
x3=(1-d)*mater_param.c3*dev(neps_cr+neps_p-neps_ii3);
x=x1+x2+x3;
sigma_eff=sigma-x;
nor_sigma_eff_dev=nor(dev(sigma_eff));
eps_cr=neps_cr+dt*(1-d)^(-mater_param.n_Norton)*mater_param.A_Norton*nor_sigma_eff_dev^(mater_param.n_Norton-1)*dev(sigma_eff);
f_visc=nor_sigma_eff_dev-sqrt(2/3)*(mater_param.K*(1-d)-mater_param.D*nf*mater_param.sigma1*exp(sigma_m/((1-nf)*mater_param.sigma1)));
lambda_p=1/mater_param.eta*max(0,f_visc);
if (lambda_p>0)
    eps_p=neps_p+dt*lambda_p*dev(sigma_eff)/nor_sigma_eff_dev;
else
    eps_p=neps_p;
end
nor_eps_i_dot=nor((eps_cr-neps_cr)/dt+(eps_p-neps_p)/dt); %lambda_i=lambda_cr+lambda_p
eps_ii1=neps_ii1+dt*mater_param.kappa1/(1-d)*nor_eps_i_dot*x1;
eps_ii2=neps_ii2+dt*mater_param.kappa2/(1-d)*nor_eps_i_dot*x2;
eps_ii3=neps_ii3+dt*mater_param.kappa3/(1-d)*nor_eps_i_dot*x3;
f=nf+dt*(1-nf)*lambda_p*mater_param.D*nf*exp(sigma_m/((1-nf)*mater_param.sigma1))+dt*mater_param.A_nuc*lambda_p;
omega=nomega+dt*mater_param.B*(1-d)^(-mater_param.m)*nor_sigma_eff_dev^mater_param.m;
end

