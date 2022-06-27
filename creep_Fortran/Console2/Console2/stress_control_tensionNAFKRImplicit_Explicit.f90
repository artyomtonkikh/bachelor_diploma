program stress_control_tensionNAFKRImplicit_Explicit
    real*8 :: k, mu, A, n_Norton, c, kappa, B, m !параметры материала
    real*8 :: sigma0 !начальное напряжение
    real*8 :: stress_current !текущее напряжение    
    real*8 :: dt, t0, tN !шаг и интервал времени
    real*8 :: nomega, omega !параметр повреждаемости damage
    real*8 :: eps(3,3), eps_cr(3,3), eps_ii(3,3), neps_cr(3,3), neps_ii(3,3) !тензоры деформаций, ползучести и ii
    real*8 :: sigma(3,3), x(3,3) !тензоры напряжений и микронапряжений
    real*8 :: Jacob(2,2) !матрица Якоби производных resid по eps
    real*8 :: inv_Jacob(2,2) !обратная матрица Jacob^(-1) к матрице Якоби
    real*8 :: deter_Jacob !определитель матрицы Якоби
    real*8 :: eps_prob1(3,3), eps_prob2(3,3) !пробные тензоры деформаций
    real*8 :: vector_eps(2) !вектор [eps11, eps22]
    real*8 :: increment(2), resid(2)
    real*8, allocatable, dimension(:) :: strain_cr(:) !массив компонент (1,1) тензора ползучести для графика
    real*8, allocatable, dimension(:) :: strain_ii(:) !массив компонент (1,1) тензора деформаций ii для графика
    real*8, allocatable, dimension(:) :: stress(:) !массив компонент (1,1) тензора напряжений для графика
    real*8, allocatable, dimension(:) :: damage(:) !массив повреждаемости omega для графика
    real*8, allocatable, dimension(:) :: t(:) !массив моментов времени для графика
    integer :: N !число шагов по времени
    integer :: n_time, z !индекс шага по времени
    integer :: alive !параметр жизни материала
    integer :: nit !число итераций для Ньютона-Рафсона
    k=80820.63d0
    mu=58525.28d0 !параметры Ламе
    A=9.3275d-15
    n_Norton=2.8801d0 !параметры закона Нортона
    c=17880.931d0
    kappa=0.01189d0
    B=2.3036d-19 !параметры из уравнения на повреждаемость (damage)
    m=5.268d0 !w'=B*||dev(sigma_eff)||^(m)*(1-w)^(-m)
    sigma0=244.0d0
    tN=4.0d6 !для 244 tN=4.0d6, для 266 tN=2.0d6
    open (777,file="parameters_set.txt",status='old')
    read (UNIT=777,FMT=*) A, n_Norton, c, kappa, B, m, sigma0, tN, N
    close (777)
    allocate (strain_cr(N+1))
    allocate (strain_ii(N+1))
    allocate (stress(N+1))
    allocate (damage(N+1))
    allocate (t(N+1))
    t0=0.0d0
    dt=(tN-t0)/N
    do n_time=1,N+1
        t(n_time)=t0+(n_time-1)*dt
    end do
    eps=0.0d0 !начальное значение
    eps_cr=0.0d0 !начальное значение
    eps_ii=0.0d0 !начальное значение
    neps_cr=0.0d0 !начальное значение
    neps_ii=0.0d0 !начальное значение
    nomega=0.0d0 !начальное значение повреждаемости
    strain_cr=0.0d0 !-->забиты 0 чтобы не было проблем<--
    strain_ii=0.0d0 !-->забиты 0 чтобы не было проблем<--
    stress=0.0d0 !-->забиты 0 чтобы не было проблем<--
    stressx=0.0d0 !-->забиты 0 чтобы не было проблем<--
    damage=0.0d0 !-->забиты 0 чтобы не было проблем<--
    alive=1
    n_time=0
    do while (n_time<(N+1) .and. alive==1)
        n_time=n_time+1
        if (n_time==181) then
            continue
        end if
        damage(n_time)=nomega
        stress_current=sigma0/(1.0d0+eps(2,2))**2.0d0
        strain_cr(n_time)=neps_cr(1,1)
        strain_ii(n_time)=neps_ii(1,1)
        call OneStepNAFKRImplicit_Explicit(eps_cr, eps_ii, sigma, x, omega, eps, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt)
        
        !метод Ньютона-Рафсона
        vector_eps(1)=eps(1,1)
        vector_eps(2)=eps(2,2)
        increment(1)=1.0d0
        increment(2)=0.0d0
        resid=0.0d0
        nit=0
        do while (sqrt(increment(1)**2.0d0+increment(2)**2.0d0)>10.0d-13 .and. nit<20)
            nit=nit+1
            Jacob=0.0d0
            resid(1)=sigma(1,1)-stress_current
            resid(2)=sigma(2,2)
            eps_prob1=0.0d0
            eps_prob2=0.0d0
            eps_prob1(1,1)=vector_eps(1)+10.0d-8 !возмущение
            eps_prob1(2,2)=vector_eps(2)
            eps_prob1(3,3)= eps_prob1(2,2)
            call OneStepNAFKRImplicit_Explicit(eps_cr, eps_ii, sigma, x, omega, eps_prob1, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt)
            Jacob(1,1)=(sigma(1,1)-stress_current-resid(1))/10.0d-8
            Jacob(2,1)=(sigma(2,2)-resid(2))/10.0d-8
            eps_prob2(1,1)=vector_eps(1)
            eps_prob2(2,2)=vector_eps(2)+10.0d-8 !возмущение
            eps_prob2(3,3)=eps_prob2(2,2)
            call OneStepNAFKRImplicit_Explicit(eps_cr, eps_ii, sigma, x, omega, eps_prob2, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt)       
            Jacob(1,2)=(sigma(1,1)-stress_current-resid(1))/10.0d-8
            Jacob(2,2)=(sigma(2,2)-resid(2))/10.0d-8
            
            !обращение матрицы Якоби и умножение вектора resid на Jacob^(-1)
            deter_Jacob=Jacob(1,1)*Jacob(2,2)-Jacob(1,2)*Jacob(2,1)
            inv_Jacob(1,1)=Jacob(2,2)/deter_Jacob
            inv_Jacob(1,2)=-Jacob(1,2)/deter_Jacob
            inv_Jacob(2,1)=-Jacob(2,1)/deter_Jacob
            inv_Jacob(2,2)=Jacob(1,1)/deter_Jacob
            !матрица Jacob^(-1) найдена
            
            increment(1)=-(inv_Jacob(1,1)*resid(1)+inv_Jacob(1,2)*resid(2))
            increment(2)=-(inv_Jacob(2,1)*resid(1)+inv_Jacob(2,2)*resid(2))
            
            vector_eps=vector_eps+increment
            eps_prob2(1,1)=vector_eps(1)
            eps_prob2(2,2)=vector_eps(2)
            eps_prob2(3,3)=eps_prob2(2,2)
            call OneStepNAFKRImplicit_Explicit(eps_cr, eps_ii, sigma, x, omega, eps_prob2, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt)
            if (nit==10) then
                nit=nit
                increment=0.0d0
                do z=(n_time+1),(N+1)
                    strain_cr(z)=strain_cr(n_time)
                alive=0
                end do
            end if
        end do
        !конец метода Ньютона-Рафсона
    
        !strain_cr(n_time)=eps_cr(1,1)
        !strain_ii(n_time)=eps_ii(1,1)
        if (omega>0.84999d0) then
            alive=0
            do z=(n_time+1),(N+1)
                strain_cr(z)=strain_cr(n_time)    
            end do
        end if
        eps(1,1)=vector_eps(1)
        eps(2,2)=vector_eps(2)
        eps(3,3)=eps(2,2)
        stress(n_time)=sigma(1,1)
        neps_cr=eps_cr !для следующего шага по времени
        neps_ii=eps_ii
        nomega=omega
    end do
    do n_time=1,N+1
        write(UNIT=777,FMT=*) t(n_time), " ", strain_cr(n_time)
    end do
end program stress_control_tensionNAFKRImplicit_Explicit