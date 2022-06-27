subroutine OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt)
    !input: eps, neps_cr, neps_p, neps_ii1, neps_ii2, nomega, nf, k, mu, c1, c2, kappa1, kappa2, n_Norton, A_Norton, K, D, sigma1, eta, A_nuc, B, m, dt
    !output: eps_cr, eps_p, eps_i, eps_ii1, eps_ii2, sigma, f, omega
    real*8 :: xi, xi_per, F_par, F_par_per, xi_left, xi_right
    real*8 :: increment, F_par_left, F_par_right, F_par_ini
    real*8 :: residuum, residuum_per, residuum_left, residuum_right
    real*8 :: derive
    real*8 :: k, mu, c1, c2, kappa1, kappa2, n_Norton, A_Norton, n_x, K_yield, D, sigma1, eta, A_nuc, B, m, dt, eta1, eta2
    real*8 :: k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, A_Norton_init, K_yield_init
    real*8 :: eps(3,3), neps_i(3,3), neps_ii1(3,3), neps_ii2(3,3)
    real*8 :: eps_i(3,3), eps_ii1(3,3), eps_ii2(3,3), sigma(3,3)
    real*8 :: eps_i_left(3,3), eps_ii1_left(3,3), eps_ii2_left(3,3), eps_i_right(3,3), eps_ii1_right(3,3), eps_ii2_right(3,3)
    real*8 :: eye(3,3), dev_eps(3,3), dev_eps_e(3,3), dev_eps_ie1(3,3), dev_eps_ie2(3,3), x(3,3), x1(3,3), x2(3,3)
    real*8 :: sigma_eff(3,3), dev_sigma_eff(3,3)
    real*8 :: damage, sigma_m, f_visc, lambda_p, lambda_cr, f, omega, nf, nomega
    real*8 :: nor, nor_sigma_eff_dev, trace, resid_function
    integer :: nit
    nit=0
    eye=0.0d0
    do i=1,3
        eye(i,i)=1.0d0
    end do
    damage=sqrt(nomega**2+nf**2)
    k=k_init*(1-damage)
    mu=mu_init*(1-damage)
    A_Norton=A_Norton_init*(1-damage)**(-n_Norton)
    c1=c1_init*(1-damage)
    c2=c2_init*(1-damage)
    kappa1=kappa1_init/(1.0d0-damage)
    kappa2=kappa2_init/(1.0d0-damage)
    call dev(dev_eps_e,eps-neps_i) !вычисление девиатора упругой деформации
    sigma=k*trace(eps-neps_i)*eye+2.0d0*mu*dev_eps_e
    sigma_m=1.0d0/3.0d0*trace(sigma)
    K_yield=K_yield_init*(1.0d0-damage)-D*nf*sigma1*dexp(sigma_m/((1.0d0-nf)*sigma1))
    call dev(dev_eps_ie1,neps_i-neps_ii1)
    call dev(dev_eps_ie2,neps_i-neps_ii2)
    x1=c1*dev_eps_ie1
    x2=c2*dev_eps_ie2
    x=x1+x2
    sigma_eff=sigma-x
    call dev(dev_sigma_eff,sigma_eff)
    nor_sigma_eff_dev=nor(dev_sigma_eff)
    F_par_ini=nor_sigma_eff_dev !начальное приближение для F
    F_par=F_par_ini
    increment=10.0d0 !инкремент F равен 10, чтобы запустить цикл метода Ньютона
    xi=dt*A_Norton*F_par**n_Norton+dt/eta*max(0.0d0,F_par-sqrt(2.0d0/3.0d0)*K_yield)
    call dev(dev_eps, eps)
    F_par=max(F_par,1.0d-5)
    F_par_left=1.0d-5
    F_par_right=F_par_ini
    residuum=resid_function(F_par, xi, dev_eps, eps_i, eps_ii1, eps_ii2, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    residuum=10.0d0
    do while (abs(residuum)>1.0d-5 .and. nit<30)
        nit=nit+1
        if (nit>10 .or. F_par<0.0d0) then
            !nit=max(nit,11) !единожды попав в бисекцию остаться в бисекции
            xi_left=dt*A_Norton*F_par_left**n_Norton+dt/eta*max(0.0d0,F_par_left-sqrt(2.0d0/3.0d0)*K_yield)
            residuum_left=resid_function(F_par_left, xi_left, dev_eps, eps_i_left, eps_ii1_left, eps_ii2_left, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt)
            xi_right=dt*A_Norton*F_par_right**n_Norton+dt/eta*max(0.0d0,F_par_right-sqrt(2.0d0/3.0d0)*K_yield)
            residuum_right=resid_function(F_par_right, xi_right, dev_eps, eps_i_right, eps_ii1_right, eps_ii2_right, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt)
            F_par=(F_par_left+F_par_right)/2.0d0
            xi=dt*A_Norton*F_par**n_Norton+dt/eta*max(0.0d0,F_par-sqrt(2.0d0/3.0d0)*K_yield)
            residuum=resid_function(F_par, xi, dev_eps, eps_i, eps_ii1, eps_ii2, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt)
            
            !работает
            if (residuum_right*residuum<0.0d0) then
                F_par_left=F_par
            else
                F_par_right=F_par
            end if
            
            !было и не работало
            !if (residuum_left*residuum<0.0d0) then
            !    F_par_right=F_par
            !end if
            !if (residuum_left*residuum>0.0d0) then
            !    F_par_left=F_par
            !end if
        else
            xi=dt*A_Norton*F_par**n_Norton+dt/eta*max(0.0d0,F_par-sqrt(2.0d0/3.0d0)*K_yield)
            F_par_per=F_par+1.0d-4
            xi_per=dt*A_Norton*F_par_per**n_Norton+dt/eta*max(0.0d0,F_par_per-sqrt(2.0d0/3.0d0)*K_yield)
            residuum=resid_function(F_par, xi, dev_eps, eps_i, eps_ii1, eps_ii2, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt) !базовое значение невязки
            residuum_per=resid_function(F_par_per, xi_per, dev_eps, eps_i, eps_ii1, eps_ii2, neps_i, neps_ii1, neps_ii2, k, mu, A_Norton, n_Norton, kappa1, kappa2, c1, c2, eta1, eta2, dt) !возмущенное значение невязки
            if(abs(residuum_per)<1.0d-10) then
                F_par=F_par_per
                exit
            end if
            derive=(residuum_per-residuum)/1.0d-4 !производная невязки по кси
            increment=-residuum_per/derive
            F_par=F_par+increment
        end if
    end do
    !F найдена
    F_par=max(F_par,1.0d-5) !для избежания деления на 0
    xi=dt*A_Norton*F_par**n_Norton+dt/eta*max(0.0d0,F_par-sqrt(2.0d0/3.0d0)*K_yield)
    call Subalgorithm_2(eps_i, dev_eps, neps_i, neps_ii1, neps_ii2, xi, F_par, mu, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    call Subalgorithm_1(eps_ii1, eps_ii2, eps_i, neps_ii1, neps_ii2, xi, kappa1, kappa2, c1, c2, eta1, eta2, dt)
    call dev(dev_eps_e,eps-eps_i)
    sigma=k*trace(eps-eps_i)*eye+2.0d0*mu*dev_eps_e
    sigma_m=1.0d0/3.0d0*trace(sigma)
    f_visc=F_par-sqrt(2.0d0/3.0d0)*K_yield
    lambda_p=1.0d0/eta*max(0.0d0,f_visc)
    if (lambda_p>0.0d0) then
        continue
    end if
    f=nf+dt*(1.0d0-nf)*lambda_p*D*nf*dexp(sigma_m/((1.0d0-nf)*sigma1))+dt*A_nuc*lambda_p
    omega=nomega+dt*B*(1.0d0-damage)**(-m)*F_par**m
end subroutine OneStep