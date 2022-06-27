real*8 function resid_function(ksi, eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt)
    !input: ksi, eps, neps_cr, neps_ii, k, mu, A, n, kappa, c, dt
    !output: resid_function
    real*8 :: ksi, k, mu, A, n_Norton, kappa, c, dt, F
    real*8 :: eps(3,3), neps_cr(3,3), neps_ii(3,3)
    real*8 :: sigma(3,3), x(3,3), sigma_eff(3,3)
    real*8 :: eps_cr(3,3), eps_ii(3,3)
    real*8 :: eye(3,3) !единичный тензор
    real*8 :: dev_sigma_eff(3,3), dev_eps_e(3,3), dev_eps_ie(3,3)
    real*8 :: trace_eps_e
    real*8 :: nor_dev_sigma_eff
    real*8 :: nor, trace
    integer :: i,j
    eye=0.0d0
    do i=1,3
        eye(i,i)=1.0d0
    end do
    F=(ksi/(A*dt))**(1.0d0/n_Norton)
    F=max(F,10.0d-10)
    call UpdateEpscr(eps_cr, eps, neps_cr, neps_ii, ksi, kappa, c, mu, F)
    call UpdateEpsii(eps_ii, neps_ii, eps_cr, ksi, kappa, c)
    trace_eps_e=trace(eps-eps_cr)
    call dev(dev_eps_e, eps-eps_cr)
    call dev(dev_eps_ie, eps_cr-eps_ii)
    sigma=k*trace_eps_e*eye+2.0d0*mu*dev_eps_e!актуальное напряжение
    x=c*dev_eps_ie !актуальное микронапряжение
    sigma_eff=sigma-x
    call dev(dev_sigma_eff, sigma_eff)
    nor_dev_sigma_eff=nor(dev_sigma_eff)
    resid_function=(ksi/(A*dt))-(nor_dev_sigma_eff)**n_Norton
end function resid_function