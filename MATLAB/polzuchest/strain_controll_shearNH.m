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
eps=zeros(3);
neps_cr=zeros(3);
dt=0.01;
t0=0;
tN=100;
N=(tN-t0)/dt;
t=t0:dt:tN;
normstrain=zeros(size(t));
normstress=zeros(size(t));
for n=1:1:(N+1)
    eps(1,2)=t(n)/tN*(0.03);
    eps(2,1)=eps(1,2);
    normstrain(n)=nor(neps_cr); %nor(A) = sqrt(trace(A*A'))
    [eps_cr, sigma]=OneStepNortonHooke(eps, neps_cr, k, mu, A, n_Norton, dt);
    normstress(n)=nor(sigma);
    neps_cr=eps_cr; %��� ���������� ���� �� �������
end
figure;
plot(t,normstrain,'r');
grid on;
figure;
plot(t,normstress,'b');
grid on;