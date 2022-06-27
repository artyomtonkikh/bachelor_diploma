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
A=8.174e-15;
n_Norton=2.9;
damping=2e-6;
nu=k/(2*(k+mu));
kappa=0.012;
c=16027;
eps=zeros(3);
neps_cr=zeros(3);
neps_ii=zeros(3);
dt=0.0005;
t0=0;
tN=20;
N=(tN-t0)/dt;
t=t0:dt:tN;
stress_hold=250; %напряжение в тесте
strain_cr=zeros(size(t));
strain_ii=zeros(size(t));
stress=zeros(size(t));
stressx=zeros(size(t));
stress_trans=zeros(size(t));
%norm=0;
%norm2=0;
for n=1:1:(N+1)
    stress_current=244/(1+eps(2,2))^2;
    %for i=1:1:3
    % for j=1:1:3
    %      norm=norm+neps_cr(i,j)^2;
    %   end
    %end
    %norm=sqrt(norm);
    %normstrain(n)=norm;
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    [eps_cr, eps_ii, sigma,x] = OneStepNortonArmstrongFrederickImplicit(eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
    eps(1,1)=eps(1,1)-(sigma(1,1)-stress_current)*damping;
    eps(2,2)=eps(2,2)-sigma(2,2)*damping;
    eps(3,3)=eps(2,2);
    %for i=1:1:3
    %   for j=1:1:3
    %      norm2=norm2+sigma(i,j)^2;
    % end
    %end
    %norm2=sqrt(norm2);
    stress(n)=sigma(1,1);
    stress_trans(n)=sigma(2,2);
    stressx(n)=x(1,1);
    neps_cr=eps_cr; %для следующего шага по времени
    neps_ii=eps_ii;
    %norm=0;
    %norm2=0;
end
figure;
plot(t,strain_cr,'r');
xlabel('t, hr');
ylabel('eps cr');
grid on;
figure;
plot(t,strain_ii,'g');
xlabel('t, hr');
ylabel('eps ii');
grid on;
figure;
plot(t,stress,'b');
xlabel('t, hr');
ylabel('stress, MPa');
grid on;
figure;
plot(t,stress_trans,'k');
xlabel('t, hr');
ylabel('stress trans, MPa');
grid on;