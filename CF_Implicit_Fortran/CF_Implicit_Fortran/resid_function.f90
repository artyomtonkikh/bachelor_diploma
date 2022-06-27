real*8 function resid_function(F_par, xi, dev_eps, eps_i, eps_ii1, eps_ii2, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    !input: F_par, xi, dev_eps, eps_i, eps_ii1, eps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt
    !output: resid_function
    real*8 :: F_par, xi, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt
    real*8 :: eps_i(3,3), dev_eps(3,3), eps_ii1(3,3), eps_ii2(3,3), neps_i(3,3), neps_ii1(3,3), neps_ii2(3,3)
    real*8 :: nor, nor_sigma_eff_deviator
    call Subalgorithm_2(eps_i, dev_eps, neps_i, neps_ii1, neps_ii2, xi, F_par, mu, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    call Subalgorithm_1(eps_ii1, eps_ii2, eps_i, neps_ii1, neps_ii2, xi, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    nor_sigma_eff_deviator=nor(2.0d0*mu*(dev_eps-eps_i)-c1*(eps_i-eps_ii1)-c2*(eps_i-eps_ii2))
    resid_function=F_par-nor_sigma_eff_deviator
end function resid_function