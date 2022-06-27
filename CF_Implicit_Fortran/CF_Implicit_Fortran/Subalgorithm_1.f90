subroutine Subalgorithm_1(eps_ii1, eps_ii2, eps_i, neps_ii1, neps_ii2, xi, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    !input: eps_i, neps_ii1, neps_ii2, xi, kappa1, kappa2, c1, c2, eta1, eta2, dt
    !output: eps_ii1, eps_ii2
    real*8 :: xi, kappa1, kappa2, c1, c2, eta1, eta2, dt
    real*8 :: eps_i(3,3), eps_ii1(3,3), eps_ii2(3,3), neps_ii1(3,3), neps_ii2(3,3)
    real*8 :: denominator1, denominator2
    denominator1=(1+kappa1*c1*(xi+dt/eta1))
    denominator2=(1+kappa2*c2*(xi+dt/eta2))
    eps_ii1=1.0d0/denominator1 * (neps_ii1+kappa1*c1*(xi+dt/eta1)*eps_i)
    eps_ii2=1.0d0/denominator2 * (neps_ii2+kappa2*c2*(xi+dt/eta2)*eps_i)
end subroutine Subalgorithm_1