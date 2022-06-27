subroutine OneStepNAFKRImplicit_Explicit(eps_cr, eps_ii, sigma, x, omega, eps, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt)
    !Implicit по eps   Explicit по omega
    !input: eps, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt
    !output: eps_cr, eps_ii, sigma, x, omega
    real*8 :: ksi, F
    real*8 :: increment
    real*8 :: residuum, residuum_per
    real*8 :: derive
    real*8 :: omega, nomega, phi, nphi
    real*8 :: k, mu, A, n_Norton, kappa, c, B, m, dt
    real*8 :: eye(3,3) !единичный тензор
    real*8 :: eps(3,3), neps_cr(3,3), neps_ii(3,3)
    real*8 :: eps_cr(3,3), eps_ii(3,3), sigma(3,3), x(3,3)
    real*8 :: dev_eps_ie(3,3), dev_eps_e(3,3)
    real*8 :: trace_eps_e, trace_eps_ie
    real*8 :: resid_function, trace
    integer :: i,j
    ksi=0.0d0 !начальное приближение для ksi
    increment=1.0d0 !инкремент кси равен 1, чтобы запустить цикл метода Ньютона
    eye=0.0d0
    do i=1,3
        eye(i,i)=1.0d0
    end do
    do while (abs(increment)>1.0d-13)
        residuum=resid_function(ksi, eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt) !базовое значение невязки
        residuum_per=resid_function(ksi+1.0d-8, eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt) !возмущенное значение невязки
        derive=(residuum_per-residuum)/1.0d-8 !производная невязки по кси
        increment=-residuum/derive !метод Ньютона
        ksi=ksi+increment
    end do
    !кси найдена
    F=(ksi/(A*dt))**(1.0d0/n_Norton)
    F=max(F,10.0d-10) !для избежания деления на 0
    call UpdateEpscr(eps_cr, eps, neps_cr, neps_ii, ksi, kappa, c, mu, F)
    call UpdateEpsii(eps_ii, neps_ii, eps_cr, ksi, kappa, c)
    nphi=-log(1.0d0-nomega);
    phi=nphi+dt*exp(nphi)*B*F**m;
    omega=1.0d0-exp(-phi);
    !omega=min(0.85d0, omega)
    trace_eps_e=trace(eps-eps_cr)
    call dev(dev_eps_e, eps-eps_cr)
    call dev(dev_eps_ie, eps_cr-eps_ii)
    sigma=(1.0d0-omega)*(k*trace_eps_e*eye+2.0d0*mu*dev_eps_e)
    x=(1.0d0-omega)*c*dev_eps_ie
end subroutine OneStepNAFKRImplicit_Explicit