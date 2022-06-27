function [res] = Error244(A,n_Norton, c, kappa,k,mu)
format long
eps=zeros(3);
neps_cr=zeros(3);
neps_ii=zeros(3);
dt=10000;
t0=0;
tN_244=2e6;%время эксперимента до проявления поврежденности
N=(tN_244-t0)/dt; %количество шагов по времени
t=t0:dt:tN_244;
strain_cr=zeros(size(t));
strain_ii=zeros(size(t));
stress=zeros(size(t));
stressx=zeros(size(t));
stress_trans=zeros(size(t));
for n=1:(N+1)
    stress_current=244/((1+eps(2,2))^2); %учёт сужения поперечного сечения
    %stress_current=266/((1+eps(2,2))^2);
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    [eps_cr, eps_ii, sigma,x] = OneStepNortonArmstrongFrederickImplicit(eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
    
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
        [~,~, sigma,~] = OneStepNortonArmstrongFrederickImplicit(eps_prob1, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
        Jacob(1,1)=(sigma(1,1)-stress_current-resid(1))/1e-8;
        Jacob(2,1)=(sigma(2,2)-resid(2))/1e-8;
        eps_prob2(1,1)=vector_eps(1);
        eps_prob2(2,2)=vector_eps(2)+1e-8;
        eps_prob2(3,3)=vector_eps(2)+1e-8;
        [~, ~, sigma,~] = OneStepNortonArmstrongFrederickImplicit(eps_prob2, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
        Jacob(1,2)=(sigma(1,1)-stress_current-resid(1))/1e-8;
        Jacob(2,2)=(sigma(2,2)-resid(2))/1e-8;
        increment=-Jacob\resid;
        vector_eps=vector_eps+increment;
        eps_prob2(1,1)=vector_eps(1);
        eps_prob2(2,2)=vector_eps(2);
        eps_prob2(3,3)=vector_eps(2);
        [eps_cr, eps_ii, sigma,x] = OneStepNortonArmstrongFrederickImplicit(eps_prob2, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
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
end
plot(t,strain_cr,':k');
grid on;
hold on;
N_exp=25;
eps_exp=dlmread('1123K_244MPa.txt');
eps_cr_exp=zeros(1,N_exp);
t_exp=zeros(1,N_exp);
for kk=1:1:N_exp
    t_exp(kk)=eps_exp(kk,1);
    eps_cr_exp(kk)=eps_exp(kk,2);
end
plot(t_exp,eps_cr_exp, '--b');
res=0;
for kk=1:1:N_exp
    distance=10000000000000;
    for n=1:1:(N+1)
        if abs(t(n)-t_exp(kk))<distance
            distance=abs(t(n)-t_exp(kk));
            neighbour=n;
        end
    end
    %сосед точки kk найден
    res=res+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
end
