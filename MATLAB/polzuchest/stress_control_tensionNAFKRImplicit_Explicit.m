%программа моделирования мягкого нагружения с учётом повреждений
%Norton-Armstrong-Frederick-Kachanov-Rabotnov
format long
k=80820.63;
mu=58525.28;
A=9.3275e-15;
n_Norton=2.8801;
c=17880.931;
kappa=0.01189;
B=2.3036e-19;%константы из уравнения на поврежденность (damage)
m=5.268;%w'=B*||dev(sigma_eff)||^(m)*(1-w)^(-m)
sigma0=244;
damping=2e-6;
%для 244 dt=20000
%для 266 dt=8912.5
N=200;
t0=0;
%для 244 tN=4e6
%для 266 tN=1782500
tN=4e6;%simulation time
dt=(tN-t0)/N;
t=t0:dt:tN;
eps=zeros(3);%тензор деформации в момент t=0
neps_cr=zeros(3);%тензор деформации ползучести в момент t=0
neps_ii=zeros(3);%тензор деформации ii в момент t=0
nomega=0;%повреждения в момент t=0
strain_cr=zeros(size(t));%массив eps_cr(1,1)
strain_ii=zeros(size(t));%массив eps_ii(1,1)
stress=zeros(size(t));%массив sigma(1,1) тензора напряжений
stressx=zeros(size(t));%массив x(1,1) тензора микронапряжений
stress_trans=zeros(size(t));%массив sigma(2,2) тензора напряжений
damage=zeros(size(t));%массив поврежденности omega
n=0;
alive=1;
while n<(N+1) && alive==1
    n=n+1;
    damage(n)=nomega;
    stress_current=sigma0/(1+eps(2,2))^2;
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKRImplicit_Explicit(eps, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt);
    
    %метод Ньютона-Рафсона
    vector_eps=[eps(1,1); eps(2,2)];
    increment=[1; 0];
    resid=[0; 0];
    nit=0;
    while nor(increment)>1e-13 && nit<20
        nit=nit+1;
        Jacob=zeros(2);
        resid(1)=sigma(1,1)-stress_current;
        resid(2)=sigma(2,2);
        eps_prob1=zeros(3);
        eps_prob2=zeros(3);
        eps_prob1(1,1)=vector_eps(1)+1e-8;%возмущение тензора eps для возмущенного eps_axial
        eps_prob1(2,2)=vector_eps(2);
        eps_prob1(3,3)=vector_eps(2);
        [~,~,sigma,~,~] = OneStepNAFKRImplicit_Explicit(eps_prob1, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt);
        Jacob(1,1)=(sigma(1,1)-stress_current-resid(1))/1e-8;
        Jacob(2,1)=(sigma(2,2)-resid(2))/1e-8;
        eps_prob2(1,1)=vector_eps(1);
        eps_prob2(2,2)=vector_eps(2)+1e-8;
        eps_prob2(3,3)=vector_eps(2)+1e-8;
        [~,~,sigma,~,~] = OneStepNAFKRImplicit_Explicit(eps_prob2, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt);
        Jacob(1,2)=(sigma(1,1)-stress_current-resid(1))/1e-8;
        Jacob(2,2)=(sigma(2,2)-resid(2))/1e-8;
        increment=-Jacob\resid;
        vector_eps=vector_eps+increment;
        eps_prob2(1,1)=vector_eps(1);
        eps_prob2(2,2)=vector_eps(2);
        eps_prob2(3,3)=vector_eps(2);
        [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKRImplicit_Explicit(eps_prob2, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt);
        if nit==19
            nit=nit;
            increment=0; %образец разрушился
            strain_cr(n+1:N+1)=strain_cr(n);
            alive=0;
        end
    end
    %конец метода Ньютона-Рафсона
%     strain_cr(n)=eps_cr(1,1);
%     strain_ii(n)=eps_ii(1,1);
    if omega>0.84999
        alive=0;
        strain_cr(n+1:N+1)=strain_cr(n);
    end
    eps(1,1)=vector_eps(1);
    eps(2,2)=vector_eps(2);
    eps(3,3)=eps(2,2);
    stress(n)=sigma(1,1);
    stress_trans(n)=sigma(2,2);
    stressx(n)=x(1,1);
    neps_cr=eps_cr; %для следующего шага по времени
    neps_ii=eps_ii;
    nomega=omega;
end
%figure;
plot(t,strain_cr);
hold on;
xlabel('t, hr');
ylabel('eps cr');
grid on;
eps_cr_fortran=dlmread('C:\Users\artyo\Documents\Visual Studio 2015\Projects\Console2\Console2\strain_cr_fortran.txt');
eps_cr_simul_fortran=zeros(1,N+1);
t_simul_fortran=zeros(1,N+1);
for kk=1:1:N+1
    t_simul_fortran(kk)=eps_cr_fortran(kk,1);
    eps_cr_simul_fortran(kk)=eps_cr_fortran(kk,2);
end
pause(1);
plot(t_simul_fortran, eps_cr_simul_fortran);
% figure;
% plot(t,strain_ii);
% xlabel('t, hr');
% ylabel('eps ii');
% grid on;
% figure;
% plot(t,stress);
% xlabel('t, hr');
% ylabel('stress, MPa');
% grid on;
% figure;
% plot(t,stress_trans);
% xlabel('t, hr');
% ylabel('stress trans, MPa');
% grid on;
% figure;
% plot(t,damage);
% xlabel('t, hr');
% ylabel('damage');
% grid on;