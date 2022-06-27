program simulation
    real*8 :: E_init, k_init, mu_init, A_Norton_init, n_Norton, n_x, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, B, m, K_yield_init, D, sigma1, eta, A_nuc, f_critical
    real*8 :: strain_rate, stress_creep, strain_max, damping, dt, dt_ep, dt_creep, t, t_max, dt_cr_ini
    real*8 :: nomega, nf, stress_max, stress_min, eps_cr_start, eps_cr_end, t_cr_start, t_cr_end, delta_eps, eps_old
    integer :: stadia, change_stadia, N_cycles, n, number_cycles, step_in_cycle, n_step, step_number_of_new_cycle
    real*8 :: eps_i(3,3), eps_ii1(3,3), eps_ii2(3,3), sigma(3,3), f, omega
    real*8 :: eps(3,3), neps_i(3,3), neps_ii1(3,3), neps_ii2(3,3)
    real*8 :: Jacob_Matrix(2,2), inv_Jacob_Matrix(2,2), deter_Jacob, Jacob_scalar, increment_vector(2), increment_scalar, resid_vector(2), resid_scalar
    real*8 :: eps_prob1(3,3), eps_prob2(3,3), eps_vector(2), eps_scalar
    real*8 :: time, stress, strain
    real*8, allocatable, dimension(:) :: stress_fr(:)
    real*8, allocatable, dimension(:) :: strain_fr(:)
    real*8, allocatable, dimension(:) :: stress_amplitude(:)
    real*8, allocatable, dimension(:) :: creep_rate(:)
    real*8, allocatable, dimension(:) :: cycles(:)
    real*8, allocatable, dimension(:) :: time_fr(:)
    integer :: fr, frame, nit
    !computation too long
    strain_rate=3600.0d-3
    fr=7
    !n_step=2000
    stadia=1
    change_stadia=0

    k_init=119837.44d0      !объёмный модуль упругости
    mu_init=58525.26d0      !модуль сдвига
    !E_init=9.0d0*k*mu/(3.0d0*k+mu)

    eta=1.0d0/strain_rate
    f_critical=0.02d0

    open (777,file="parameters_set.txt",status='old')
    read (UNIT=777,FMT=*) A_Norton_init, n_Norton, n_x, c1_init, c2_init, c3_init, kappa1_init, kappa2_init, kappa3_init, eta1, eta2, eta3, B, m, K_yield_init, D, sigma1, A_nuc, stress_creep, strain_max, number_cycles, n_step
    
    dt=1.0d-5/5.0d0
    dt_ep=dt
    !dt_cr_ini=4.640679d-4
    dt_cr_ini=1.0d-4
    dt_creep=dt_cr_ini
    t=0.0d0;
    t_max=0.2d0             !simulation time
    nomega=0.0d0            !повреждения при ползучести в момент t=0
    nf=0.0d0                !повреждения при пластичности в момент t=0

    eps=0.0d0
    neps_i=0.0d0
    neps_ii1=0.0d0
    neps_ii2=0.0d0
    
    n=0
    step_in_cycle=0
    stress_max=0.0d0
    stress_min=0.0d0
    eps_cr_start=0.0d0
    eps_cr_end=strain_max
    t_cr_start=0.0d0
    t_cr_end=0.0d0
    N_cycles=1
    
    open(1, file= "sigma_eps.txt")
    open(2, file= "creep_rate.txt")
    open(3, file= "stress_amplitude.txt")
    
    allocate (stress_fr(1000000))
    stress=0.0d0
    allocate (strain_fr(1000000))
    strain=0.0d0
    allocate (stress_amplitude(1000000))
    stress_amplitude=0.0d0
    allocate (creep_rate(1000000))
    creep_rate=0.0d0
    allocate (cycles(1000000))
    cycles=0.0d0
    allocate (time_fr(1000000))
    time=0.0d0
    frame=0
    step_number_of_new_cycle=0
    do while (N_cycles<=number_cycles .and. step_in_cycle<n_step)
        n=n+1
        step_in_cycle=step_in_cycle+1
        strain=eps(1,1)
        if (stadia==1) then
            if(n==1) then
               time=dt_ep
            else if (n/=1) then
                time=time+dt_ep
            end if
            eps(1,1)=eps(1,1)+strain_rate*dt_ep !elasto-plastic tension
            call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
            !метод Ньютона
            eps_scalar=eps(2,2)
            increment_scalar=1.0d0
            resid_scalar=2.0d0
            nit=0
            do while (abs(resid_scalar)>1.0d-13 .and. nit<20)
                nit=nit+1
                eps_prob2=0.0d0
                eps_prob2(1,1)=eps(1,1)
                eps_prob2(2,2)=eps_scalar
                eps_prob2(3,3)=eps_prob2(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob2, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
                resid_scalar=sigma(2,2)
                eps_prob2(1,1)=eps(1,1)
                eps_prob2(2,2)=eps_scalar+1.0d-5 !возмущение
                eps_prob2(3,3)=eps_prob2(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob2, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
                Jacob_scalar=(sigma(2,2)-resid_scalar)/1.0d-5
                increment_scalar=-resid_scalar/Jacob_scalar
                eps_scalar=eps_scalar+increment_scalar
            end do
            eps(2,2)=eps_scalar
            eps(3,3)=eps(2,2)
            call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
            if (sigma(1,1)>=stress_creep) then
                change_stadia=1
                dt_creep=dt_cr_ini
                if (N_cycles==1) then
                    eps_cr_start=eps(1,1)
                end if
                t_cr_start=time
            end if
        end if
        if (stadia==2) then
            time=time+dt_creep
            stress_current=stress_creep
            if(time>22) then
                continue
            end if
            !метод Ньютона-Рафсона
            eps_vector(1)=eps(1,1)
            eps_vector(2)=eps(2,2)
            eps_old=eps(1,1) !запомнили деформацию с прошлого шага, чтобы вычислить приращение
            increment_vector(1)=1.0d0
            increment_vector(2)=0.0d0
            resid_vector(1)=2.0d0
            resid_vector(2)=2.0d0
            nit=0
            do while (sqrt(resid_vector(1)**2.0d0+resid_vector(2)**2.0d0)>1.0d-6 .and. nit<50)
                nit=nit+1
                Jacob_Matrix=0.0d0
                eps_prob1=0.0d0
                eps_prob1(1,1)=eps_vector(1)
                eps_prob1(2,2)=eps_vector(2)
                eps_prob1(3,3)=eps_prob1(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob1, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_creep)
                resid_vector(1)=sigma(1,1)-stress_current
                resid_vector(2)=sigma(2,2)
                eps_prob1(1,1)=eps_vector(1)+1.0d-4 !возмущение
                eps_prob1(2,2)=eps_vector(2)
                eps_prob1(3,3)=eps_prob1(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob1, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_creep)
                Jacob_Matrix(1,1)=(sigma(1,1)-stress_current-resid_vector(1))/1.0d-4
                Jacob_Matrix(2,1)=(sigma(2,2)-resid_vector(2))/1.0d-4
                eps_prob2=0.0d0
                eps_prob2(1,1)=eps_vector(1)
                eps_prob2(2,2)=eps_vector(2)+1.0d-4 !возмущение
                eps_prob2(3,3)=eps_prob2(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob2, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_creep)
                Jacob_Matrix(1,2)=(sigma(1,1)-stress_current-resid_vector(1))/1.0d-4
                Jacob_Matrix(2,2)=(sigma(2,2)-resid_vector(2))/1.0d-4
            
                !обращение матрицы Якоби и умножение вектора resid на Jacob^(-1)
                deter_Jacob=Jacob_Matrix(1,1)*Jacob_Matrix(2,2)-Jacob_Matrix(1,2)*Jacob_Matrix(2,1)
                inv_Jacob_Matrix(1,1)=Jacob_Matrix(2,2)/deter_Jacob
                inv_Jacob_Matrix(1,2)=-Jacob_Matrix(1,2)/deter_Jacob
                inv_Jacob_Matrix(2,1)=-Jacob_Matrix(2,1)/deter_Jacob
                inv_Jacob_Matrix(2,2)=Jacob_Matrix(1,1)/deter_Jacob
                !матрица Jacob^(-1) найдена
            
                increment_vector(1)=-(inv_Jacob_Matrix(1,1)*resid_vector(1)+inv_Jacob_Matrix(1,2)*resid_vector(2))
                increment_vector(2)=-(inv_Jacob_Matrix(2,1)*resid_vector(1)+inv_Jacob_Matrix(2,2)*resid_vector(2))
            
                eps_vector=eps_vector+increment_vector
                if(eps_vector(1)/=eps_vector(1)) then
                    continue
                end if
            end do
            
            eps(1,1)=eps_vector(1)
            delta_eps=eps(1,1)-eps_old
            eps(2,2)=eps_vector(2)
            eps(3,3)=eps(2,2)
            call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_creep)
            if ((nit>5) .or. (delta_eps>1.0d-5)) then
                dt_creep=dt_creep*0.9d0
            end if
            if((nit<3) .and. (delta_eps<1.0d-5)) then
                dt_creep=dt_creep*1.1d0
            end if
            !dt_creep=min(1.0d0, dt_creep)
            if (eps(1,1)>=strain_max) then
              change_stadia=1
              t_cr_end=time
              eps_cr_end=eps(1,1)
              !A_Norton_init=A_Norton_init/3.0d0
              if(stress_max<sigma(1,1)) then
                    stress_max=sigma(1,1)
              end if
              creep_rate(N_cycles)=(eps_cr_end-eps_cr_start)/(t_cr_end-t_cr_start)
            end if
        end if
        if (stadia==3) then
            time=time+dt_ep
            eps(1,1)=eps(1,1)-strain_rate*dt_ep
            !метод Ньютона
            eps_scalar=eps(2,2)
            increment_scalar=1.0d0
            resid_scalar=2.0d0
            nit=0
            do while (abs(resid_scalar)>1.0d-10 .and. nit<20)
                nit=nit+1
                eps_prob2=0.0d0
                eps_prob2(1,1)=eps(1,1)
                eps_prob2(2,2)=eps_scalar
                eps_prob2(3,3)=eps_prob2(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob2, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
                resid_scalar=sigma(2,2)
                eps_prob2(1,1)=eps(1,1)
                eps_prob2(2,2)=eps_scalar+1.0d-5 !возмущение
                eps_prob2(3,3)=eps_prob2(2,2)
                call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps_prob2, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
                Jacob_scalar=(sigma(2,2)-resid_scalar)/1.0d-5
                increment_scalar=-resid_scalar/Jacob_scalar
                eps_scalar=eps_scalar+increment_scalar
            end do
            eps(2,2)=eps_scalar
            eps(3,3)=eps(2,2)
            call OneStep(eps_i, eps_ii1, eps_ii2, sigma, f, omega, eps, neps_i, neps_ii1, neps_ii2, nomega, nf, k_init, mu_init, c1_init, c2_init, kappa1_init, kappa2_init, eta1, eta2, n_Norton, A_Norton_init, n_x, K_yield_init, D, sigma1, eta, A_nuc, B, m, dt_ep)
            !if (eps(1,1)<-1.0d-4) then
            !    continue
            !end if
            if (eps(1,1)<=-strain_max) then
               change_stadia=1
                stress_min=sigma(1,1)
           end if
        end if
        stress=sigma(1,1)
        neps_i=eps_i !для следующего шага по времени
        neps_ii1=eps_ii1
        neps_ii2=eps_ii2
        neps_ii3=eps_ii3
        nomega=omega
        nf=f
        if (change_stadia==1) then
           stadia=stadia+1
           if (stadia==4) then
               stadia=1
           end if
           change_stadia=0
        end if
        if (sqrt(nomega**2+(nf/f_critical)**2)>=0.8d0) then
            exit
        end if
        !write(1,*) eps(1,1), sigma(1,1), time
        if (int(n/fr)*fr==n) then
            frame=frame+1
            stress_fr(frame)=stress
            time_fr(frame)=time
            strain_fr(frame)=strain
            write(1,*) strain_fr(frame), stress_fr(frame), time_fr(frame)
        end if
        if (eps(1,1)>=0.0d0 .and. strain<=0.0d0 .and. n>1 .and. abs(n-step_number_of_new_cycle)>=10) then
            stress_amplitude(N_cycles)=(stress_max-stress_min)/2.0d0
            N_cycles=N_cycles+1
            cycles(N_cycles-1)=N_cycles-1
            step_in_cycle=0
            step_number_of_new_cycle=n
        end if
    end do
    if (step_in_cycle==n_step) then
        write(1,*) 'computation too long'
        write(2,*) 'computation too long'
        write(3,*) 'computation too long'
    else if (sqrt(nomega**2+(nf/f_critical)**2)>=0.5d0) then
        write(1,*) 'destruction'
        write(2,*) 'destruction'
        write(3,*) 'destruction'
    else
        do i=1,frame
           !write(1,*) strain_fr(i), stress_fr(i), time_fr(i)
        end do
        do i=1,N_cycles-1
            write(2,*) cycles(i), creep_rate(i)
            write(3,*) cycles(i), stress_amplitude(i)
        end do
    end if
    close(1)
    close(2)
    close(3)
end program simulation