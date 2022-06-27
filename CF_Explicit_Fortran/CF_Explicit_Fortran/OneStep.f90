subroutine OneStep(eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega, eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, k, mu, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, n_Norton, A_Norton, n_x, K_yield, D, sigma1, eta, A_nuc, B, m, dt)
    !input: eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, k, mu, c1, c2, c3, kappa1, kappa2, kappa3, n_Norton, A_Norton, K, D, sigma1, eta, A_nuc, B, m, dt
    !output: eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega
    real*8 :: k, mu, c1, c2, c3, kappa1, kappa2, kappa3, n_Norton, A_Norton, n_x, K_yield, D, sigma1, eta, A_nuc, B, m, dt, eta1, eta2, eta3
    real*8 :: eps(3,3), neps_cr(3,3), neps_p(3,3), neps_ii1(3,3), neps_ii2(3,3), neps_ii3(3,3)
    
    real*8 :: eps_cr(3,3), eps_p(3,3), eps_ii1(3,3), eps_ii2(3,3), eps_ii3(3,3), sigma(3,3)
    
    real*8 :: eye(3,3), dev_eps_e(3,3), dev_eps_ie1(3,3), dev_eps_ie2(3,3), dev_eps_ie3(3,3), x(3,3), x1(3,3), x2(3,3), x3(3,3)
    real*8 :: sigma_eff(3,3), dev_sigma_eff(3,3)
    real*8 :: damage, sigma_m, f_visc, lambda_p, f, omega, nf, nomega
    real*8 :: nor, nor_sigma_eff_dev, nor_eps_i_dot, trace
    eye=0.0d0
    do i=1,3
        eye(i,i)=1.0d0
    end do
    damage=sqrt(nomega**2+nf**2)
    call dev(dev_eps_e,eps-neps_cr-neps_p)
    sigma=(1.0d0-damage)*(k*trace(eps-neps_cr-neps_p)*eye+2.0d0*mu*dev_eps_e)
    sigma_m=1.0d0/3.0d0*trace(sigma)
    call dev(dev_eps_ie1,neps_cr+neps_p-neps_ii1)
    call dev(dev_eps_ie2,neps_cr+neps_p-neps_ii2)
    call dev(dev_eps_ie3,neps_cr+neps_p-neps_ii3)
    x1=(1.0d0-damage)*c1*dev_eps_ie1
    x2=(1.0d0-damage)*c2*dev_eps_ie2
    x3=(1.0d0-damage)*c3*dev_eps_ie3
    x=x1+x2+x3
    sigma_eff=sigma-x
    call dev(dev_sigma_eff,sigma_eff)
    nor_sigma_eff_dev=nor(dev_sigma_eff)
    eps_cr=neps_cr+dt*(1.0d0-damage)**(-n_Norton)*A_Norton*nor_sigma_eff_dev**(n_Norton-1.0d0)*dev_sigma_eff
    f_visc=nor_sigma_eff_dev-sqrt(2.0d0/3.0d0)*(K_yield*(1.0d0-damage)-D*nf*sigma1*exp(sigma_m/((1.0d0-nf)*sigma1)))
    lambda_p=1.0d0/eta*max(0.0d0,f_visc)
    if (lambda_p>0) then
        eps_p=neps_p+dt*lambda_p*dev_sigma_eff/nor_sigma_eff_dev
    else if (lambda_p<=0) then
        eps_p=neps_p
    end if
    nor_eps_i_dot=nor((eps_cr-neps_cr)/dt+(eps_p-neps_p)/dt) !lambda_i=lambda_cr+lambda_p
    nor_x1=nor(x1)
    nor_x2=nor(x2)
    nor_x3=nor(x3)
    eps_ii1=neps_ii1+dt*kappa1/(1.0d0-damage)*(nor_eps_i_dot+(1.0d0/eta1)*nor_x1**(n_x-1.0d0))*x1
    eps_ii2=neps_ii2+dt*kappa2/(1.0d0-damage)*(nor_eps_i_dot+(1.0d0/eta2)*nor_x2**(n_x-1.0d0))*x2
    eps_ii3=neps_ii3+dt*kappa3/(1.0d0-damage)*(nor_eps_i_dot+(1.0d0/eta3)*nor_x3**(n_x-1.0d0))*x3
    f=nf+dt*(1.0d0-nf)*lambda_p*D*nf*exp(sigma_m/((1.0d0-nf)*sigma1))+dt*A_nuc*lambda_p
    omega=nomega+dt*B*(1.0d0-damage)**(-m)*nor_sigma_eff_dev**m
end subroutine OneStep