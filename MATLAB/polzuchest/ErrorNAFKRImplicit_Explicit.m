function [res] = ErrorNAFKRImplicit_Explicit(y)
format long
A=9.3275e-15*exp(y(1)-1);
n_Norton=2.8801*y(2);
c=17880.931*exp(y(3)-1);
kappa=0.01189*exp(y(4)-1);
B=2.3036e-19*exp(y(5)-1);
m=5.268*y(6);
k=80820.63;
mu=58525.28;
%рассматриваем начальное напряжение 244 MPa
sigma0=244;
t0=0;
tN244=4e6;
N=1000;
dt=(tN244-t0)/N;
t=t0:dt:tN244;
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
%     strain_cr(n)=neps_cr(1,1);
%     strain_ii(n)=neps_ii(1,1);
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
    while nor(increment)>1e-12 && nit<20
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
            strain_cr(n+1:N+1)=strain_cr(n);%strain_cr(n:N+1)=strain_cr(n-1);
            alive=0;
        end
    end
    %конец метода Ньютона-Рафсона
    
    %strain_cr(n)=eps_cr(1,1);
    %strain_ii(n)=eps_ii(1,1);
    if omega>0.84999
        alive=0;
        strain_cr(n+1:N+1)=strain_cr(n);%strain_cr(n:N+1)=strain_cr(n-1);
    end
%     if omega>0.85
%         strain_cr(n+1:N+1)=1.5*strain_cr(n);
%         break;
%     end
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
plot(t,strain_cr,':k');
hold on;
N_exp=50;
eps_exp=dlmread('1123K_244MPa.txt');
eps_cr_exp=zeros(1,N_exp);
t_exp=zeros(1,N_exp);
for kk=1:1:N_exp
    t_exp(kk)=eps_exp(kk,1);
    eps_cr_exp(kk)=eps_exp(kk,2);
end
plot(t_exp,eps_cr_exp, '--b');
res1=0;
for kk=1:1:N_exp
    distance=10000000000000;
    for n=1:1:(N+1)
        if abs(t(n)-t_exp(kk))<distance
            distance=abs(t(n)-t_exp(kk));
            neighbour=n;
        end
    end
    %сосед точки kk найден
    res1=res1+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
%рассматриваем начальное напряжение 266 MPa
sigma0=266;
tN266=2e6;
N=1000;
dt=(tN266-t0)/N;
t=t0:dt:tN266;
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
    while nor(increment)>1e-15 && nit<20
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
            strain_cr(n+1:N+1)=strain_cr(n);%strain_cr(n:N+1)=strain_cr(n-1);
            alive=0;
        end
    end
    %конец метода Ньютона-Рафсона
%     strain_cr(n)=eps_cr(1,1);
%     strain_ii(n)=eps_ii(1,1);
    if omega>0.84999
        alive=0;
        strain_cr(n+1:N+1)=strain_cr(n);%strain_cr(n:N+1)=strain_cr(n-1);
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
plot(t,strain_cr,'-.r');
N_exp=50;
eps_exp=dlmread('1123K_266MPa.txt');
eps_cr_exp=zeros(1,N_exp);
t_exp=zeros(1,N_exp);
for kk=1:1:N_exp
    t_exp(kk)=eps_exp(kk,1);
    eps_cr_exp(kk)=eps_exp(kk,2);
end
plot(t_exp,eps_cr_exp,'-g');
grid on;
xlabel('t, h');
ylabel('\epsilon_{cr}');
hold off;
pause(0.1);
res2=0;
for kk=1:1:N_exp
    distance=10000000000000;
    for n=1:1:(N+1)
        if abs(t(n)-t_exp(kk))<distance
            distance=abs(t(n)-t_exp(kk));
            neighbour=n;
        end
    end
    %сосед точки kk найден
    res2=res2+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
legend({'num 244','exp 244','num 266','exp 266'},'Location','northwest');
res=res1+res2;
end

