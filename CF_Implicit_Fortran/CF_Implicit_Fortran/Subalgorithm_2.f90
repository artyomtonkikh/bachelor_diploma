subroutine Subalgorithm_2(eps_i, dev_eps, neps_i, neps_ii1, neps_ii2, xi, F_par, mu, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    !input: dev_eps, neps_i, neps_ii1, neps_ii2, xi, F, mu, kappa1, kappa2, c1, c2, eta1, eta2, dt
    !output: eps_i
    real*8 :: xi, F_par, mu, kappa1, kappa2, c1, c2, eta1, eta2, dt
    real*8 :: eps_i(3,3), dev_eps(3,3), neps_i(3,3), neps_ii1(3,3), neps_ii2(3,3)
    real*8 :: denominator, b1, b2
    real*8 :: A1(3,3), A2(3,3)
    A1=neps_ii1/(1.0d0+kappa1*c1*(xi+dt/eta1))
    A2=neps_ii2/(1.0d0+kappa2*c2*(xi+dt/eta2))
    b1=kappa1*c1*(xi+dt/eta1)/(1.0d0+kappa1*c1*(xi+dt/eta1))
    b2=kappa2*c2*(xi+dt/eta2)/(1.0d0+kappa2*c2*(xi+dt/eta2))
    denominator=1.0d0+xi/F_par*(2.0d0*mu+c1+c2-c1*b1-c2*b2)
    eps_i=1.0d0/denominator * (neps_i+xi/F_par*(2.0d0*mu*dev_eps+c1*A1+c2*A2))
end subroutine Subalgorithm_2