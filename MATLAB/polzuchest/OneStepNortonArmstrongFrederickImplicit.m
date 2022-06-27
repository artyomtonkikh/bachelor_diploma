function [eps_cr, eps_ii, sigma, x] = OneStepNortonArmstrongFrederickImplicit(eps, neps_cr, neps_ii, k, mu, A, n, kappa, c, dt)
ksi=0;%начальное приближение для ksi
increment=1;%инкремент кси равен 1, чтобы запустить цикл метода Ньютона
while abs(increment)>1e-13
   residuum=resid(ksi,eps, neps_cr, neps_ii,  k, mu, A, n, kappa, c, dt);%базовое значение невязки
   residuum_per=resid(ksi+1e-8,eps, neps_cr, neps_ii,  k, mu, A, n, kappa, c, dt);%возмущенное значение невязки
   derive=(residuum_per-residuum)/1e-8;%производная невяки по кси
   increment=-residuum/derive;%метод Ньютона
   ksi=ksi+increment;
end
%кси найдена
F=(ksi/(A*dt))^(1/n);
F=max(F,1e-10);%для избежания деления на 0
eps_cr=UpdateEpscr(eps, neps_cr, neps_ii, ksi,kappa, c, mu, F);
eps_ii=UpdateEpsii(neps_ii, eps_cr,ksi,kappa,c);
sigma=k*trace(eps-eps_cr)*eye(3)+2*mu*dev(eps-eps_cr);%актуальное напряжение
x=c*dev(eps_cr-eps_ii); %актуальное микронапряжение
end