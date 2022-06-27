function [eps_cr, sigma] = OneStepNortonHooke(eps, neps_cr, k, mu, A, n, dt)
%������� ��� ���������� ������ ���� �� ������� ��� ������ �������-����
%input: eps-���������� ������ ����������
%neps_cr - ������ ������ ���������� ����������
%k, mu, A, n - ��������� ���������

sigma=k*trace(eps-neps_cr)*eye(3)+2*mu*dev(eps-neps_cr);
norm_stress_dev=nor(dev(sigma));
eps_cr_dot=norm_stress_dev^(n-1)*dev(sigma);%�������� ���������� ����������
eps_cr=neps_cr+dt*A*eps_cr_dot;%������������ ������� ���������� ����������
sigma=k*trace(eps-eps_cr)*eye(3)+2*mu*dev(eps-eps_cr);
end