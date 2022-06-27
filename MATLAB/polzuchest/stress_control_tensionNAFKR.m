%��������� ������������� ������� ���������� � ������ �����������
%Norton-Armstrong-Frederick-Kachanov-Rabotnov
format long
k0=80820.63;
mu0=58525.28;
A0=8.4399e-15;
n_Norton=2.8801;
c0=16609;
kappa0=0.0116;
B=2.3036e-19;%��������� �� ��������� �� �������������� (damage)
m=5.268;%w'=B*||dev(sigma_eff)||^(m)*(1-w)^(-m)
sigma0=244;
damping=2e-6;
%��� 244 dt=1000
%��� 266 dt=445.625
dt=125;
t0=0;
%��� 244 tN=4e6
%��� 266 tN=1782500
tN=4e6;%simulation time
N=(tN-t0)/dt;
t=t0:dt:tN;
eps=zeros(3);%������ ���������� � ������ t=0
neps_cr=zeros(3);%������ ���������� ���������� � ������ t=0
neps_ii=zeros(3);%������ ���������� ii � ������ t=0
nomega=0;%����������� � ������ t=0
strain_cr=zeros(size(t));%������ eps_cr(1,1)
strain_ii=zeros(size(t));%������ eps_ii(1,1)
stress=zeros(size(t));%������ sigma(1,1) ������� ����������
stressx=zeros(size(t));%������ x(1,1) ������� ���������������
stress_trans=zeros(size(t));%������ sigma(2,2) ������� ����������
damage=zeros(size(t));%������ �������������� omega
for n=1:(N+1)
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    damage(n)=nomega;
    stress_current=sigma0/(1+eps(2,2))^2;
    [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKR(eps, neps_cr, neps_ii, nomega, k0, mu0, A0, n_Norton, kappa0, c0, B, m, dt);
    if omega>0.85
        strain_cr(n+1:N+1)=strain_cr(n);
        break;
    end
    eps(1,1)=eps(1,1)-(sigma(1,1)-stress_current)*damping;
    eps(2,2)=eps(2,2)-sigma(2,2)*damping;
    eps(3,3)=eps(2,2);
    stress(n)=sigma(1,1);
    stress_trans(n)=sigma(2,2);
    stressx(n)=x(1,1);
    neps_cr=eps_cr; %��� ���������� ���� �� �������
    neps_ii=eps_ii;
    nomega=omega;
end
plot(t,strain_cr);
xlabel('t, hr');
ylabel('eps cr');
grid on;
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