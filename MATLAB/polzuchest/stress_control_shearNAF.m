%{
параметры A и n модели Нортона определены из документа stevens1981.pdf
из таблицы на странице 7, строки 244 MPa и 266 MPa (все при 1123 К) без разрушения

параметры Ламе закона Гука определены
из документа in_738alloy_preliminarydata_497.pdf
через значения модуля Юнга E и коэффициента Пуассона v
из таблицы на странице 4, строка  1600 °F
%}
format long
k=80820.63;
mu=58525.28;
A=7e-18;
n_Norton=6;
damping=1e-6;
nu=k/(2*(k+mu));
kappa=1/1000;
c=80000;
eps=zeros(3);
neps_cr=zeros(3);
neps_ii=zeros(3);
dt=0.001;
t0=0;
tN=10;
N=(tN-t0)/dt;
t=t0:dt:tN;
stress_hold=250; %напряжение в тесте
normstrain_cr=zeros(size(t));
normstrain_ii=zeros(size(t));
normstress=zeros(size(t));
stressx=zeros(size(t));
for n=1:1:(N+1)
    if mod(fix(t(n)),2)==0
        stress_current=300;
    end
    if mod(fix(t(n)),2)==1
        stress_current=100;
    end
    normstrain_cr(n)=neps_cr(1,2); %nor(A) = sqrt(trace(A*A'))
    normstrain_ii(n)=neps_ii(1,2);
    [eps_cr, eps_ii, sigma,x] = OneStepNortonArmstrongFrederick(eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
    if sigma(1,2)-x(1,2)<-50
        k=k;
    end
    eps(1,2)=eps(1,2)-(sigma(1,2)-stress_current)*damping;
    eps(2,1)=eps(1,2);
    normstress(n)=sigma(1,2);
    stressx(n)=x(1,2);
    neps_cr=eps_cr; %для следующего шага по времени
    neps_ii=eps_ii;
end
figure;
plot(t,normstrain_cr,'r');
xlabel('t, hr');
ylabel('eps cr');
grid on;
figure;
plot(t,normstrain_ii,'g');
xlabel('t, hr');
ylabel('eps ii');
grid on;
figure;
plot(t,normstress,'b');
xlabel('t, hr');
ylabel('stress, MPa');
grid on;
figure;
plot(t,stressx,'k');
xlabel('t, hr');
ylabel('x');
grid on;