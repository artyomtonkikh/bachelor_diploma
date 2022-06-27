function [res] = Error_CF_amps_RFD(mater_param)
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
strain_max=[0.38756625776162346/100 0.5902014235953355/100 0.973345449038316/100];
number_cycles=3;
n_step=100000;
res=[];
load('CF_amps.txt');
strain_amplitude_exp=CF_amps(:,1);
stress_amplitude_exp=CF_amps(:,2);
stress_amplitude_num=[];

% strain amp = 0.40016766248018243%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(1), number_cycles, n_step);
!CF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(1)=stress_amplitude(number_cycles,2);
    res(1)=stress_amplitude_exp(1)-stress_amplitude_num(1);
end


% strain amp = 0.5000152420436507%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(2), number_cycles, n_step);
!CF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(2)=stress_amplitude(number_cycles,2);
    res(2)=stress_amplitude_exp(2)-stress_amplitude_num(2);
end


% strain amp = 0.5997988050237754%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep, strain_max(3), number_cycles, n_step);
!CF_Explicit_Fortran.exe
amp=fopen('stress_amplitude.txt');
str=fgetl(amp);
fclose(amp);
if(str==" computation too long" || str==" destruction" || str==" computation too long" || str==" destruction")
    res=100000000000;
    return;
else
    load('stress_amplitude.txt');
    stress_amplitude_num(3)=stress_amplitude(number_cycles,2);
    res(3)=stress_amplitude_exp(3)-stress_amplitude_num(3);
end
%res=10*res;
figure(3);
plot(strain_amplitude_exp,stress_amplitude_exp,'-b', strain_amplitude_exp, stress_amplitude_num, '--dr');
grid on;
xlabel('strain amplitude, %');
ylabel('stress amplirude, MPa');
legend({'exp','num'},'Location','northwest');
pause(0.001);
end

