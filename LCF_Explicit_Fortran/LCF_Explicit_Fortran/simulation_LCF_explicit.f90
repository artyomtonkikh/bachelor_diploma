program simulation
    real*8 :: E, k, mu, A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, eta, A_nuc, f_critical
    real*8 :: strain_rate, stress_creep, strain_max, damping, dt, dt_ep, dt_creep, t, t_max
    real*8 :: nomega, nf, stress_max, stress_min, eps_cr_start, eps_cr_end, t_cr_start, t_cr_end
    integer :: stadia, change_stadia, N_cycles, n, number_cycles, step_in_cycle, n_step, step_number_of_new_cycle
    real*8 :: eps_cr(3,3), eps_p(3,3), eps_ii1(3,3), eps_ii2(3,3), eps_ii3(3,3), sigma(3,3), f, omega
    real*8 :: eps(3,3), neps_cr(3,3), neps_p(3,3), neps_ii1(3,3), neps_ii2(3,3), neps_ii3(3,3)
    real*8 :: time, stress, strain
    real*8, allocatable, dimension(:) :: stress_fr(:)
    real*8, allocatable, dimension(:) :: strain_fr(:)
    real*8, allocatable, dimension(:) :: stress_amplitude(:)
    real*8, allocatable, dimension(:) :: creep_rate(:)
    real*8, allocatable, dimension(:) :: cycles(:)
    real*8, allocatable, dimension(:) :: time_fr(:)
    integer :: fr, frame
    !computation too long
    strain_rate=3600.0d-3
    fr=10
    !n_step=2000
    stadia=1
    change_stadia=0

    k=119837.44d0      !объёмный модуль упругости
    mu=58525.26d0      !модуль сдвига
    E=9.0d0*k*mu/(3.0d0*k+mu)

    !A_Norton=16.0d0*3600.0d0*8.4399d-13 !константы закона Нортона
    !n_Norton=2.8801d0

    !c1=16609.0d0
    !c2=11000.0d0/2.0d0
    !c3=12000.0d0
    !kappa1=0.0116d0
    !kappa2=0.01d0
    !kappa3=0.02d0

    !B=200.0d0*3600.0d0*2.3036d-19  !константы из уравнения
    !m=5.268d0                  !повреждаемости при ползучести

    !K_yield=350.0d0         !начальный предел текучести
    !D=2.0d-2                !константы из уравнения повреждаемости при пластичности
    !sigma1=100.0d0
    eta=1.0d0/strain_rate
    !A_nuc=1.0d-5
    f_critical=0.02d0

    open (777,file="parameters_set.txt",status='old')
    read (UNIT=777,FMT=*) A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max, number_cycles, n_step
    damping=1.0d0/E*0.2d0         !? надо будет пристреливаться
    dt=1.0d-5/2.0d0
    dt_ep=dt
    dt_creep=dt
    t=0.0d0;
    t_max=0.2d0             !simulation time
    nomega=0.0d0            !повреждения при ползучести в момент t=0
    nf=0.0d0                !повреждения при пластичности в момент t=0

    eps=0.0d0
    neps_cr=0.0d0
    neps_p=0.0d0
    neps_ii1=0.0d0
    neps_ii2=0.0d0
    neps_ii3=0.0d0
    
    n=0
    step_in_cycle=0
    stress_max=stress_creep
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
        if (eps(1,1)>=0.0d0 .and. strain<=0.0d0 .and. n>1 .and. abs(n-step_number_of_new_cycle)>=10) then
            N_cycles=N_cycles+1
            cycles(N_cycles-1)=N_cycles-1
            step_in_cycle=0
            step_number_of_new_cycle=n
        end if
        n=n+1
        step_in_cycle=step_in_cycle+1
        strain=eps(1,1)
        if (stadia==1) then
            if(n==1) then
                time=dt_ep
            else if (n/=1) then
                time=time+dt_ep
            end if
            call OneStep(eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega, eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, k, mu, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, n_Norton, A_Norton, n_x, K_yield, D, sigma1, eta, A_nuc, B, m, dt_ep)
            eps(1,1)=eps(1,1)+strain_rate*dt_ep !elasto-plastic tension
            eps(2,2)=eps(2,2)-sigma(2,2)*damping
            eps(3,3)=eps(2,2)
            if (eps(1,1)>=strain_max) then
                change_stadia=1
                stress_max=sigma(1,1)
            end if
        end if
        if (stadia==2) then
            time=time+dt_ep
            call OneStep(eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega, eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, k, mu, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3,  n_Norton, A_Norton, n_x, K_yield, D, sigma1, eta, A_nuc, B, m, dt_ep)
            eps(1,1)=eps(1,1)-strain_rate*dt_ep
            eps(2,2)=eps(2,2)-sigma(2,2)*damping
            eps(3,3)=eps(2,2)
            if (eps(1,1)<=-strain_max) then
               change_stadia=1
               stress_min=sigma(1,1)
               stress_amplitude(N_cycles)=(stress_max-stress_min)/2.0d0
           end if
        end if
        stress=sigma(1,1)
        neps_cr=eps_cr !для следующего шага по времени
        neps_p=eps_p
        neps_ii1=eps_ii1
        neps_ii2=eps_ii2
        neps_ii3=eps_ii3
        nomega=omega
        nf=f
        if (change_stadia==1) then
           stadia=stadia+1
           if (stadia==3) then
               stadia=1
               !write(*,*) n
           end if
           change_stadia=0
        end if
        if (sqrt(nomega**2+(nf/f_critical)**2)>=0.8d0) then
            exit
        end if
        if (int(n/fr)*fr==n) then
            frame=frame+1
            stress_fr(frame)=stress
            time_fr(frame)=time
            strain_fr(frame)=strain
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
            write(1,*) strain_fr(i), stress_fr(i), time_fr(i)
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