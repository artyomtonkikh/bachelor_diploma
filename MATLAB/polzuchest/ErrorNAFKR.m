function [res] = ErrorNAFKR(y)
format long
A0=8.17e-15;
n_Norton=2.9;
c0=16027;
kappa0=0.0124;
B=y(1)*2.248e-19;
m=y(2)*5.289;
k0=80820.63;
mu0=58525.28;
%������������� ��������� ���������� 244 MPa
sigma0=244;
damping=2e-6;
dt=1000;
t0=0;
tN244=4e6;
N=(tN244-t0)/dt;
t=t0:dt:tN244;
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
    stress_current=sigma0/(1+eps(2,2))^2;
    [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKR(eps, neps_cr, neps_ii, nomega, k0, mu0, A0, n_Norton, kappa0, c0, B, m, dt);
    if omega>0.85
        strain_cr(n+1:N+1)=strain_cr(n);
        break;
    end
    eps(1,1)=eps(1,1)-(sigma(1,1)-stress_current)*damping;
    eps(2,2)=eps(2,2)-sigma(2,2)*damping;
    eps(3,3)=eps(2,2);
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    damage(n)=nomega;
    stress(n)=sigma(1,1);
    stress_trans(n)=sigma(2,2);
    stressx(n)=x(1,1);
    neps_cr=eps_cr; %��� ���������� ���� �� �������
    neps_ii=eps_ii;
    if omega>=1
        nomega=1;
    else
        nomega=omega;
    end
end
plot(t,strain_cr,':k');
hold on;
grid on;
N_exp=50;
eps_exp=dlmread('1123K_244MPa.txt');
eps_cr_exp=zeros(1,N_exp);
t_exp=zeros(1,N_exp);
for kk=1:1:N_exp
    t_exp(kk)=eps_exp(kk,1);
    eps_cr_exp(kk)=eps_exp(kk,2);
end
plot(t_exp,eps_cr_exp,'--b');
res1=0;
for kk=1:1:N_exp
    distance=10000000000000;
    for n=1:1:(N+1)
        if abs(t(n)-t_exp(kk))<distance
            distance=abs(t(n)-t_exp(kk));
            neighbour=n;
        end
    end
    %����� ����� kk ������
    res1=res1+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
%������������� ��������� ���������� 266 MPa
sigma0=266;
dt=445.625;
tN266=1782500;
N=(tN266-t0)/dt;
t=t0:dt:tN266;
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
    stress_current=sigma0/(1+eps(2,2))^2;
    [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKR(eps, neps_cr, neps_ii, nomega, k0, mu0, A0, n_Norton, kappa0, c0, B, m, dt);
    if omega>0.85
        strain_cr(n+1:N+1)=strain_cr(n);
        break;
    end
    eps(1,1)=eps(1,1)-(sigma(1,1)-stress_current)*damping;
    eps(2,2)=eps(2,2)-sigma(2,2)*damping;
    eps(3,3)=eps(2,2);
    strain_cr(n)=neps_cr(1,1);
    strain_ii(n)=neps_ii(1,1);
    damage(n)=nomega;
    stress(n)=sigma(1,1);
    stress_trans(n)=sigma(2,2);
    stressx(n)=x(1,1);
    neps_cr=eps_cr; %��� ���������� ���� �� �������
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
    %����� ����� kk ������
    res2=res2+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
res=1.5*res1+res2;
legend({'num 244','exp 244','num 266','exp 266'},'Location','northwest');
res=res+0;
end