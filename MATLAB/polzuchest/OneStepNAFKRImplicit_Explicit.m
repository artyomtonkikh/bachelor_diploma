function [eps_cr, eps_ii, sigma, x, omega] = OneStepNAFKRImplicit_Explicit(eps, neps_cr, neps_ii, nomega, k0, mu0, A0, n_Norton, kappa0, c0, B, m, dt)
%Implicit �� eps
%Explicit �� omega
ksi=0;%��������� ����������� ��� ksi
increment=1;%��������� ��� ����� 1, ����� ��������� ���� ������ �������
while abs(increment)>1e-12
   residuum=resid(ksi, eps, neps_cr, neps_ii, k0, mu0, A0, n_Norton, kappa0, c0, dt);%������� �������� �������
   residuum_per=resid(ksi+1e-8, eps, neps_cr, neps_ii, k0, mu0, A0, n_Norton, kappa0, c0, dt);%����������� �������� �������
   derive=(residuum_per-residuum)/1e-8;%����������� ������ �� ���
   increment=-residuum/derive;%����� �������
   ksi=ksi+increment;
end
%��� �������
F=(ksi/(A0*dt))^(1/n_Norton);
F=max(F,1e-10);%��� ��������� ������� �� 0
eps_cr=UpdateEpscr(eps, neps_cr, neps_ii, ksi, kappa0, c0, mu0, F);
eps_ii=UpdateEpsii(neps_ii, eps_cr, ksi, kappa0, c0);
nphi=-log(1-nomega);
phi=nphi+dt*exp(nphi)*B*F^m;
omega=1-exp(-phi);
%omega=min(0.85, omega); %�����������, ��� damage<=1
sigma=(1-omega)*(k0*trace(eps-eps_cr)*eye(3)+2*mu0*dev(eps-eps_cr));%���������� ����������
x=(1-omega)*c0*dev(eps_cr-eps_ii); %���������� ���������������
end

