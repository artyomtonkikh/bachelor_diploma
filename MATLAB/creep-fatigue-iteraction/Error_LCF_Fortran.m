function [res] = Error_LCF_Fortran(mater_param)
A_Norton=mater_param.A_Norton;
n_Norton=mater_param.n_Norton;
n_x=mater_param.n_x;
c1=mater_param.c1;
c2=mater_param.c2;
c3=mater_param.c3;
kappa1=mater_param.kappa1;
kappa2=mater_param.kappa2;
kappa3=mater_param.kappa3;
eta1=mater_param.eta1;
eta2=mater_param.eta2;
eta3=mater_param.eta3;
B=mater_param.B;
m=mater_param.m;
K_yield=mater_param.K_yield;
D=mater_param.D;
sigma1=mater_param.sigma1;
A_nuc=mater_param.A_nuc;
stress_creep=230;
strain_max=[0.40016766248018243/100 0.5000152420436507/100 0.5997988050237754/100 0.8002286306547965/100 1.0005365199365919/100];
number_cycles=15;
n_step=100000;
res=0;
load('sigma_eps_amplitude.txt');
strain_amplitude_exp=sigma_eps_amplitude(:,1);
stress_amplitude_exp=sigma_eps_amplitude(:,2);
stress_amplitude_num=[];

% strain amp = 0.40016766248018243%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(1), number_cycles, n_step);
!LCF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(1)=stress_amplitude(number_cycles,2);
    res=res+(stress_amplitude_exp(1)-stress_amplitude_num(1))^2;
end


% strain amp = 0.5000152420436507%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(2), number_cycles, n_step);
!LCF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(2)=stress_amplitude(number_cycles,2);
    res=res+(stress_amplitude_exp(2)-stress_amplitude_num(2))^2;
end


% strain amp = 0.5997988050237754%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(3), number_cycles, n_step);
!LCF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(3)=stress_amplitude(number_cycles,2);
    res=res+(stress_amplitude_exp(3)-stress_amplitude_num(3))^2;
end


% strain amp = 0.8002286306547965%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(4), number_cycles, n_step);
!LCF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(4)=stress_amplitude(number_cycles,2);
    res=res+(stress_amplitude_exp(4)-stress_amplitude_num(4))^2;
end


% strain amp = 1.0005365199365919%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(5), number_cycles, n_step);
!LCF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(5)=stress_amplitude(number_cycles,2);
    res=res+(stress_amplitude_exp(5)-stress_amplitude_num(5))^2;
end


figure(3);
plot(strain_amplitude_exp,stress_amplitude_exp,'b', strain_amplitude_exp, stress_amplitude_num, '-sr');
grid on;
xlabel('strain amplitude, %');
ylabel('stress amplirude, MPa');
legend({'exp','num', 'cf exp'},'Location','northwest');
pause(0.1)
%res=10*res;
end

