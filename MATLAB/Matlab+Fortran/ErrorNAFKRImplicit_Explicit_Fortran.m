function [res] = ErrorNAFKRImplicit_Explicit_Fortran(y)
format long
A=9.3275e-15*exp(y(1)-1); %9.3275e-15
n_Norton=2.88008*y(2); %2.8801
c=18797*exp(y(3)-1); %17880.931
kappa=0.0115*exp(y(4)-1); %0.01189
B=2.30358e-19*exp(y(5)-1); %2.3036e-19
m=5.27*y(6); %5.268
sigma0=244;
tN=4e6;
N=200;
change_parameters(A, n_Norton, c, kappa, B, m, sigma0, tN, N);
!Fortran_NAFKR.exe
load('fort.777');
strain_cr = fort(:,2);
t=fort(:,1);
plot(t,strain_cr,'k');
hold on;
N_exp=50;
eps_exp=dlmread('1123K_244MPa.txt');
eps_cr_exp=zeros(1,N_exp);
t_exp=zeros(1,N_exp);
for kk=1:1:N_exp
    t_exp(kk)=eps_exp(kk,1);
    eps_cr_exp(kk)=eps_exp(kk,2);
end
plot(t_exp,eps_cr_exp);
res1=0;
for kk=1:1:N_exp
    distance=10000000000000;
    for n=1:1:(N+1)
        if abs(t(n)-t_exp(kk))<distance
            distance=abs(t(n)-t_exp(kk));
            neighbour=n;
        end
    end
    %сосед точки kk найден
    res1=res1+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
%рассматриваем начальное напряжение 266 MPa
sigma0=266;
tN=2e6;
change_parameters(A, n_Norton, c, kappa, B, m, sigma0, tN, N);
!Fortran_NAFKR.exe
load('fort.777');
strain_cr = fort(:,2);
t=fort(:,1);
plot(t,strain_cr,'k');
N_exp=50;
eps_exp=dlmread('1123K_266MPa.txt');
eps_cr_exp=zeros(1,N_exp);
t_exp=zeros(1,N_exp);
for kk=1:1:N_exp
    t_exp(kk)=eps_exp(kk,1);
    eps_cr_exp(kk)=eps_exp(kk,2);
end
plot(t_exp,eps_cr_exp);
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
    %сосед точки kk найден
    res2=res2+(eps_cr_exp(kk)-strain_cr(neighbour))^2; 
end
res=res1+res2;
end

