%{
��������� A � n ������ ������� ���������� �� ��������� stevens1981.pdf
�� ������� �� �������� 7, ������ 244 MPa � 266 MPa (��� ��� 1123 �) ��� ����������

��������� ���� ������ ���� ����������
�� ��������� in_738alloy_preliminarydata_497.pdf
����� �������� ������ ���� E � ������������ �������� v
�� ������� �� �������� 4, ������  1600 �F
%}
format long
k=80820.63;
mu=58525.28;
A=(7.23e-23*3600*1e+12)^2;
n_Norton=5.88931;
nu=k/(2*(k+mu));
kappa=1e-6;
c=30000;
eps=zeros(3);
neps_cr=zeros(3);
neps_ii=zeros(3);
dt=0.01;
t0=0;
tN=100;
N=(tN-t0)/dt;
t=t0:dt:tN;
normstrain_cr=zeros(size(t));
normstrain_ii=zeros(size(t));
normstress=zeros(size(t));
for n=1:1:(N+1)
    eps(1,2)=t(n)/tN*(0.03);
    eps(2,1)=eps(1,2);
    normstrain_cr(n)=nor(neps_cr); %nor(A) = sqrt(trace(A*A'))
    normstrain_ii(n)=nor(neps_ii);
    [eps_cr, eps_ii, sigma] = OneStepNortonArmstrongFrederick(eps, neps_cr, neps_ii, k, mu, A, n_Norton, kappa, c, dt);
    normstress(n)=nor(sigma);
    neps_cr=eps_cr; %��� ���������� ���� �� �������
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