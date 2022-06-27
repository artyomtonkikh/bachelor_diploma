%программа для моделирования мягкого нагружения NH
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
A=6.7746e-18;
n_Norton=5.88931;
damping=1e-6;
nu=k/(2*(k+mu));
eps=zeros(3);
neps_cr=zeros(3);
dt=0.001;
t0=0;
tN=10;
N=(tN-t0)/dt;
t=t0:dt:tN;
stress_hold=250; %напряжение в тесте
normstrain=zeros(size(t));
normstress=zeros(size(t));
for n=1:1:(N+1)
    if mod(fix(t(n)),2)==0
        stress_current=300;
    end
    if mod(fix(t(n)),2)==1
        stress_current=200;
    end
    normstrain(n)=nor(neps_cr); %nor(A) = sqrt(trace(A*A'))
    [eps_cr, sigma]=OneStepNortonHooke(eps, neps_cr, k, mu, A, n_Norton, dt);
    eps(1,2)=eps(1,2)-(sigma(1,2)-stress_current)*damping;
    eps(2,1)=eps(1,2);
    normstress(n)=sigma(1,2);
    neps_cr=eps_cr; %для следующего шага по времени
end
figure;
plot(t,normstrain,'r');
grid on;
figure;
plot(t,normstress,'b');
grid on;