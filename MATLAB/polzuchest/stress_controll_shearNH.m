%��������� ��� ������������� ������� ���������� NH
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
stress_hold=250; %���������� � �����
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
    neps_cr=eps_cr; %��� ���������� ���� �� �������
end
figure;
plot(t,normstrain,'r');
grid on;
figure;
plot(t,normstress,'b');
grid on;