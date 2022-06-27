%программа для моделирования мягкого нагружения NAF Implicit
%{
параметры A и n модели Нортона определены из документа stevens1981.pdf
из таблицы на странице 7, строки 244 MPa и 266 MPa (все при 1123 К) без
прерывания

параметры Ламе закона Гука определены
из документа in_738alloy_preliminarydata_497.pdf
через значения модуля Юнга E и коэффициента Пуассона v
из таблицы на странице 4, строка  1600 °F
%}
%A 1e-12 n 2 c 63977 kappa 7e-3
format long
k=80820.63;
mu=58525.28;
A=8.174e-15;
n_Norton=2.9;
damping=2e-6;
nu=k/(2*(k+mu));
kappa=0.0124;
c=16027;
eps=zeros(3);
neps_cr=zeros(3);
neps_ii=zeros(3);
dt=50000;
t0=0;
tN=1000000;
N=(tN-t0)/dt;
t=t0:dt:tN;
strain_cr=zeros(size(t));
strain_ii=zeros(size(t));
stress=zeros(size(t));
stressx=zeros(size(t));
stress_trans=zeros(size(t));
for n=1:1:(N+1)
    stress_current=244/((1+eps(2,2))^2);
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

figure;
plot(t,strain_cr,'r');
xlabel('t, s');
ylabel('eps cr num');
grid on;
figure;
plot(t,strain_ii,'g');
xlabel('t, s');
ylabel('eps ii');
grid on;
figure;
plot(t,stress,'b');
xlabel('t, s');
ylabel('stress, MPa');
grid on;
figure;
plot(t,stress_trans,'k');
xlabel('t, s');
ylabel('stress trans, MPa');
grid on;