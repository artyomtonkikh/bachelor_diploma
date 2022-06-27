%программа моделирования мягкого нагружения с учётом повреждений
clear;
format long
strain_rate=1e-3*3600; %strain_rate per hour
stress_creep=230;
strain_max=0.01;
stadia=1;
change_stadia=0;

mater_param.k=80820.63;         %константы Ламе
mater_param.mu=58525.28;

mater_param.A_Norton=16*8.4399e-13*3600;%константы закона Нортона
mater_param.n_Norton=2.8801;

mater_param.c1=16609;
mater_param.c2=11000/2;
mater_param.c3=12000;
mater_param.kappa1=0.0116;
mater_param.kappa2=0.01;
mater_param.kappa3=0.02;

mater_param.B=200*2.3036e-19*3600;  %константы из уравнения
mater_param.m=5.268;       %повреждаемости при ползучести

mater_param.K=350;            %начальный предел текучести %константы из уравнения повреждаемости
mater_param.D=2e-2;           %при пластичности
mater_param.sigma1=100;
mater_param.eta=1/strain_rate;
mater_param.A_nuc=1e-5;
f_critical=0.04;

damping=1e-1; %? надо будет пристреливаться
dt=5e-7/10;
dt_ep=dt;
dt_creep=dt*160;
t=0;
t_max=1;    %simulation time
nomega=0;            %повреждения при ползучести в момент t=0
nf=0;                %повреждения при пластичности в момент t=0
eps=zeros(3);
neps_cr=zeros(3);
neps_p=zeros(3);
neps_ii1=zeros(3);
neps_ii2=zeros(3);
neps_ii3=zeros(3);

n=0;
stress_max=stress_creep;
stress_min=0;
eps_cr_start=0;
eps_cr_end=strain_max;
t_cr_start=0;
t_cr_end=0;
N_cycles=0;
while(N_cycles<=1) 
    n=n+1;
    if(n>1 && eps(1,1)>=0 && strain(n-1)<=0)
       N_cycles=N_cycles+1; 
       if(N_cycles>1)
          stress_amplitude(N_cycles-1)=(stress_max-stress_min)/2;
          creep_rate(N_cycles-1)=(eps_cr_end-eps_cr_start)/(t_cr_end-t_cr_start);
          cycles(N_cycles-1)=N_cycles-1;
       end
    end
    strain(n)=eps(1,1);
%     strain_cr(n)=neps_cr(1,1);
%     strain_p(n)=neps_p(1,1);
%     damage(n)=sqrt(nomega^2+(nf/f_critical)^2);
    
%     [stress_current]=current_stress_table(t(n));
%     [eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega] = OneStep(eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, mater_param, dt);
%     eps(1,1)=eps(1,1)-(sigma(1,1)-stress_current)*damping;
%     eps(2,2)=eps(2,2)-sigma(2,2)*damping;
%     eps(3,3)=eps(2,2);
    if (stadia==1)
        if(n==1)
           t(1)=dt_ep; 
        else
            t(n)=t(n-1)+dt_ep;
        end
        [eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega] = OneStep(eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, mater_param, dt_ep);
        eps(1,1)=eps(1,1)+strain_rate*dt_ep; %elasto-plastic tension
        eps(2,2)=eps(2,2)-sigma(2,2)*damping*dt_ep;
        eps(3,3)=eps(2,2);
        if (sigma(1,1)>=stress_creep)
           change_stadia=1;
           eps_cr_start=eps(1,1);
           t_cr_start=t(n);
       end
    end
    if(stadia==2)
        t(n)=t(n-1)+dt_creep;
        stress_current=stress_creep;
        [eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega] = OneStep(eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, mater_param, dt_creep);
        eps(1,1)=eps(1,1)-(sigma(1,1)-stress_current)*damping*dt_creep;
        eps(2,2)=eps(2,2)-sigma(2,2)*damping*dt_creep;
        eps(3,3)=eps(2,2);
        if(eps(1,1)>=strain_max)
           change_stadia=1;
           t_cr_end=t(n);
        end
    end
    if(stadia==3)
        t(n)=t(n-1)+dt_ep;
        [eps_cr, eps_p, eps_ii1, eps_ii2, eps_ii3, sigma, f, omega] = OneStep(eps, neps_cr, neps_p, neps_ii1, neps_ii2, neps_ii3, nomega, nf, mater_param, dt_ep);
        eps(1,1)=eps(1,1)-strain_rate*dt_ep; 
        eps(2,2)=eps(2,2)-sigma(2,2)*damping*dt_ep;
        eps(3,3)=eps(2,2);
        if(eps(1,1)<=-strain_max)
           change_stadia=1;
           stress_min=sigma(1,1);
       end
    end
    stress(n)=sigma(1,1);
%     stress_trans(n)=sigma(2,2);
    neps_cr=eps_cr; %для следующего шага по времени
    neps_p=eps_p;
    neps_ii1=eps_ii1;
    neps_ii2=eps_ii2;
    neps_ii3=eps_ii3;
    nomega=omega;
    nf=f;
    if(change_stadia==1)
       stadia=stadia+1;
       if(stadia==4)
           stadia=1;
       end
       change_stadia=0;
    end
    if (sqrt(nomega^2+(nf/f_critical)^2)>=0.5)
        n=n;
       break; 
    end
end
%figure;
plot(strain,stress);
xlabel('eps');
ylabel('stress, MPa');
grid on;
% figure;
% semilogx(cycles,stress_amplitude);
% xlabel('N cycles');
% ylabel('stress amplitude, MPa');
% grid on;
% figure;
% loglog(cycles,creep_rate);
% xlabel('N cycles');
% ylabel('creep rate, 1/h');
% grid on;
% figure;
% plot(t,strain);
% xlabel('t');
% ylabel('eps');
% grid on;
% figure;
% plot(t,strain_cr);
% xlabel('t');
% ylabel('eps cr');
% grid on;
% figure;
% plot(t,strain_p);
% xlabel('t');
% ylabel('eps p');
% grid on;
% plot(t,stress);
% xlabel('t');
% ylabel('stress, MPa');
% grid on;
% figure;
% plot(t,damage);
% xlabel('t');
% ylabel('damage');
% grid on;