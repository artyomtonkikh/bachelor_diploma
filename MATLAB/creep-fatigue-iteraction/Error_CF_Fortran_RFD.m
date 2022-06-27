function [res] = Error_CF_Fortran_RFD(mater_param)
%clc;
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
stress_creep=[150 170 190 230];
strain_max=[0.004 0.006 0.01];
number_cycles=[10 7 2];
n_step=1000000;
res_strain1=[];
res_strain2=[];
res_strain3=[];
res_strain4=[];
res_stress1=[];
res_stress2=[];
res_stress3=[];
bad=0;
%разные деформации
% eps=0.4%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(4), strain_max(1), number_cycles(1), n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=number_cycles(1);
    for n=1:N
        res_stress1(n)=1000000;
    end
    bad=1;
else
    load('stress_amplitude.txt');
    cycles_eps_0_4=stress_amplitude(:,1);
    stress_amplitude_num_eps_0_4=stress_amplitude(:,2);
    load('stress_amplitude_0_4.txt');
    cycles_exp_eps_0_4=stress_amplitude_0_4(:,1);
    stress_amplitude_exp_eps_0_4=stress_amplitude_0_4(:,2);
    N_data=length(cycles_eps_0_4);
    for i=1:N_data
        u_eps_0_4(i)=linear_interpolation(cycles_eps_0_4(i), cycles_exp_eps_0_4, stress_amplitude_exp_eps_0_4);
    end
    N=number_cycles(1);
    for n=1:N
        data=u_eps_0_4(n);
        resid=stress_amplitude_num_eps_0_4(n)-data;
        res_stress1(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];
% eps=0.6%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(4), strain_max(2), number_cycles(2), n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=number_cycles(2);
    for n=1:N
        res_stress2(n)=1000000;
    end
    bad=1;
else
    load('stress_amplitude.txt');
    cycles_eps_0_6=stress_amplitude(:,1);
    stress_amplitude_num_eps_0_6=stress_amplitude(:,2);
    load('stress_amplitude_0_6.txt');
    cycles_exp_eps_0_6=stress_amplitude_0_6(:,1);
    stress_amplitude_exp_eps_0_6=stress_amplitude_0_6(:,2);
    N_data=length(cycles_eps_0_6);
    for i=1:N_data
        u_eps_0_6(i)=linear_interpolation(cycles_eps_0_6(i), cycles_exp_eps_0_6, stress_amplitude_exp_eps_0_6);
    end
    N=number_cycles(2);
    for n=1:N
        data=u_eps_0_6(n);
        resid=stress_amplitude_num_eps_0_6(n)-data;
        res_stress2(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];
% eps=1%
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(4), strain_max(3), number_cycles(3), n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=number_cycles(3);
    for n=1:N
        res_stress3(n)=1000000;
    end
    bad=1;
else
    load('stress_amplitude.txt');
    cycles_eps_1=stress_amplitude(:,1);
    stress_amplitude_num_eps_1=stress_amplitude(:,2);
    load('stress_amplitude_1.txt');
    cycles_exp_eps_1=stress_amplitude_1(:,1);
    stress_amplitude_exp_eps_1=stress_amplitude_1(:,2);
    N_data=length(cycles_eps_1);
    for i=1:N_data
        u_eps_1(i)=linear_interpolation(cycles_eps_1(i), cycles_exp_eps_1, stress_amplitude_exp_eps_1);
    end
    N=number_cycles(3);
    for n=1:N
        data=u_eps_1(n);
        resid=stress_amplitude_num_eps_1(n)-data;
        res_stress3(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];

%разные напряжения
%stress=150 MPa
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(1), strain_max(1), 15, n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=15;
    for n=1:N
        res_strain1(n)=1000000;
    end
    bad=1;
else
    load('creep_rate.txt');
    cycles_150_MPa=creep_rate(:,1);
    creep_rate_num_150_MPa=creep_rate(:,2);
    load('creep_rate_150_MPa.txt');
    cycles_exp_150_MPa=creep_rate_150_MPa(:,1);
    creep_rate_exp_150_MPa=creep_rate_150_MPa(:,2);
    N_data=length(cycles_150_MPa);
    for i=1:N_data
        u_150_MPa(i)=linear_interpolation(cycles_150_MPa(i), cycles_exp_150_MPa, creep_rate_exp_150_MPa);
    end
    N=15;
    for n=1:N
        data=u_150_MPa(n);
        resid=log10(creep_rate_num_150_MPa(n)/3600)-log10(data);
        res_strain1(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];
%stress=170 MPa
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(2), strain_max(1), 15, n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=15;
    for n=1:N
        res_strain2(n)=1000000;
    end
    bad=1;
else
    load('creep_rate.txt');
    cycles_170_MPa=creep_rate(:,1);
    creep_rate_num_170_MPa=creep_rate(:,2);
    load('creep_rate_170_MPa.txt');
    cycles_exp_170_MPa=creep_rate_170_MPa(:,1);
    creep_rate_exp_170_MPa=creep_rate_170_MPa(:,2);
    N_data=length(cycles_170_MPa);
    for i=1:N_data
        u_170_MPa(i)=linear_interpolation(cycles_170_MPa(i), cycles_exp_170_MPa, creep_rate_exp_170_MPa);
    end
    N=15;
    for n=1:N
        data=u_170_MPa(n);
        resid=log10(creep_rate_num_170_MPa(n)/3600)-log10(data);
        res_strain2(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];
%stress=190 MPa
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(3), strain_max(1), 15, n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=15;
    for n=1:N
        res_strain3(n)=1000000;
    end
    bad=1;
else
    load('creep_rate.txt');
    cycles_190_MPa=creep_rate(:,1);
    creep_rate_num_190_MPa=creep_rate(:,2);
    load('creep_rate_190_MPa.txt');
    cycles_exp_190_MPa=creep_rate_190_MPa(:,1);
    creep_rate_exp_190_MPa=creep_rate_190_MPa(:,2);
    N_data=length(cycles_190_MPa);
    for i=1:N_data
        u_190_MPa(i)=linear_interpolation(cycles_190_MPa(i), cycles_exp_190_MPa, creep_rate_exp_190_MPa);
    end
    N=15;
    for n=1:N
        data=u_190_MPa(n);
        resid=log10(creep_rate_num_190_MPa(n)/3600)-log10(data);
        res_strain3(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];
%stress=230 MPa
change_parameters(A_Norton, n_Norton, n_x, c1, c2, c3, kappa1, kappa2, kappa3, eta1, eta2, eta3, B, m, K_yield, D, sigma1, A_nuc, stress_creep(4), strain_max(1), 15, n_step);
!CF_Implicit_Fortran.exe
cr=fopen('creep_rate.txt');
sa=fopen('stress_amplitude.txt');
str_cr=fgetl(cr);
str_sa=fgetl(sa);
fclose(cr);
fclose(sa);
if(str_cr==" computation too long" || str_cr==" destruction" || str_sa==" computation too long" || str_sa==" destruction")
    N=15;
    for n=1:N
        res_strain4(n)=1000000;
    end
    bad=1;
else
    load('creep_rate.txt');
    cycles_230_MPa=creep_rate(:,1);
    creep_rate_num_230_MPa=creep_rate(:,2);
    load('creep_rate_230_MPa.txt');
    cycles_exp_230_MPa=creep_rate_230_MPa(:,1);
    creep_rate_exp_230_MPa=creep_rate_230_MPa(:,2);
    N_data=length(cycles_230_MPa);
    for i=1:N_data
        u_230_MPa(i)=linear_interpolation(cycles_230_MPa(i), cycles_exp_230_MPa, creep_rate_exp_230_MPa);
    end
    N=15;
    for n=1:N
        data=u_230_MPa(n);
        resid=log10(creep_rate_num_230_MPa(n)/3600)-log10(data);
        res_strain4(n)=resid;
    end
end
stress_amplitude=[];
creep_rate=[];
if(bad==1)
    a=0;
else
    figure(1);
    semilogx(cycles_eps_0_4,stress_amplitude_num_eps_0_4,'--or',cycles_eps_0_4,u_eps_0_4,'--r',cycles_eps_0_6,stress_amplitude_num_eps_0_6,'-.sb',cycles_eps_0_6,u_eps_0_6,'-.b',cycles_eps_1,stress_amplitude_num_eps_1,':^k',cycles_eps_1,u_eps_1,':k');
    xlabel('cycles');
    ylabel('stress amplitude, MPa');
    ylim([0 500])
    legend({'num \epsilon_{max}=0.4%','exp \epsilon_{max}=0.4%','num \epsilon_{max}=0.6%','exp \epsilon_{max}=0.6%','num \epsilon_{max}=1%','exp \epsilon_{max}=1%'},'Location','southwest');
    % hold on;
    % semilogx(cycles_eps_0_4,u_eps_0_4,'-ok');
    % pause(0.5);
    % 
    % semilogx(cycles_eps_0_6,stress_amplitude_num_eps_0_6,'-s');
    % hold on;
    % semilogx(cycles_eps_0_6,u_eps_0_6,'-sk');
    % pause(0.5);
    % 
    % semilogx(cycles_eps_1,stress_amplitude_num_eps_1,'-^');
    % hold on;
    % semilogx(cycles_eps_1,u_eps_1,'-^k');
    figure(2);
    loglog(cycles_150_MPa,creep_rate_num_150_MPa/3600,'-db',cycles_150_MPa,u_150_MPa,'-b',cycles_170_MPa,creep_rate_num_170_MPa/3600,'--^r',cycles_170_MPa,u_170_MPa,'--r',cycles_190_MPa,creep_rate_num_190_MPa/3600,'-.sg',cycles_190_MPa,u_190_MPa,'-.g',cycles_230_MPa,creep_rate_num_230_MPa/3600,':ok',cycles_230_MPa,u_230_MPa,':k');
    xlabel('cycles');
    ylabel('creep rate, 1/h');
    ylim([1e-10 1e-3]);
    legend({'num \sigma_h=150 MPa','exp \sigma_h=150 MPa','num \sigma_h=170 MPa','exp \sigma_h=170 MPa','num \sigma_h=190 MPa','exp \sigma_h=190 MPa','num \sigma_h=230 MPa','exp \sigma_h=230 MPa'},'Location','southeast');
    % hold on;
    % loglog(cycles_150_MPa,u_150_MPa,'-dk');
    % pause(0.5);
    % 
    % loglog(cycles_170_MPa,creep_rate_num_170_MPa/3600,'-^');
    % hold on;
    % loglog(cycles_170_MPa,u_170_MPa,'-^k');
    % pause(0.5);
    % 
    % loglog(cycles_190_MPa,creep_rate_num_190_MPa/3600,'-s');
    % hold on;
    % loglog(cycles_190_MPa,u_190_MPa,'-sk');
    % pause(0.5);
    % 
    % loglog(cycles_230_MPa,creep_rate_num_230_MPa/3600,'-o');
    % hold on;
    % loglog(cycles_230_MPa,u_230_MPa,'-ok');
end
u_eps_0_4=[];
u_eps_0_6=[];
u_eps_1=[];
u_150_MPa=[];
u_170_MPa=[];
u_190_MPa=[];
u_230_MPa=[];
%cla;
res=[res_stress1, res_stress2, res_stress3, 100*res_strain1, 100*res_strain2, 100*res_strain3, 100*res_strain4];
%res=100*(res_strain1+res_strain2+res_strain3+res_strain4)+res_stress1+1.5*res_stress2+10*res_stress3;
end

