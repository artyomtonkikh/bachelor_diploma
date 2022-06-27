clear;
clc;
format long
k=80820.63;
mu=58525.28;
A=8.17e-15;
n_Norton=2.9;
c=16027;
kappa=0.0124;
B=2.248e-19;%константы из уравнения на поврежденность (damage)
m=5.289;%w'=B*||dev(sigma_eff)||^(m)*(1-w)^(-m)
sigma0=244;
damping=2e-6;
%для 244 dt=1000
%для 266 dt=445.625
dt=2000;
t0=0;
%для 244 tN=4e6
%для 266 tN=1782500
tN=4e6;%simulation time
N=(tN-t0)/dt;
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
for n=1:(N+1)
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    damage(n)=nomega;
    stress_current=sigma0/(1+eps(2,2))^2;
    [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKRImplicit_Explicit(eps, neps_cr, neps_ii, nomega, k, mu, A, n_Norton, kappa, c, B, m, dt);
    if omega>0.85
        strain_cr(n+1:N+1)=strain_cr(n);
        break;
    end
    
    %метод Ньютона-Рафсона
    vector_eps=[eps(1,1); eps(2,2)];
    increment=[1; 0];
    resid=[0; 0];
    while nor(increment)>1e-15
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
    end
    %конец метода Ньютона-Рафсона
    
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
plot(t,strain_cr);
xlabel('t, hr');
ylabel('eps cr');
grid on;
hold on;