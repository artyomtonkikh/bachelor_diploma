function [res]=Error_3_experiments(y)
%disp(y);
res=0;
res_creep_rate=0;
res_stress=0;
res1=0;
res2=0;
%идеальные первые и вторые графики явной схемы
% mater_param.A_Norton=4.9017e-6*100^(y(1)-1);
% mater_param.n_Norton=2.5125*10^(y(2)-1);
% mater_param.eta1=1.6746e+03*100^(y(3)-1);
% mater_param.eta2=1.3561e+03*100^(y(4)-1);
% mater_param.c1=9.6884e+03*10^(y(5)-1);
% mater_param.c2=4.1459e+04*10^(y(6)-1);
% mater_param.kappa1=0.0090*exp(y(7)-1);
% mater_param.kappa2=0.0028*exp(y(8)-1);
% mater_param.K_yield=393.0147*10^(y(9)-1);
% mater_param.n_x=1;

%параметры неявной схемы
mater_param.A_Norton=5.2862e-6;%*exp(y(1)-1);
mater_param.n_Norton=2.5189;%*exp(y(2)-1);
mater_param.eta1=1.8438e3;%*exp(y(3)-1);
mater_param.eta2=1.0967e3;%*exp(y(4)-1);
mater_param.c1=1.0229e4;%*exp(y(5)-1);
mater_param.c2=4.2454e4;%*exp(y(6)-1);
mater_param.kappa1=0.0089;%*exp(y(7)-1);
mater_param.kappa2=0.0027;%*exp(y(8)-1);
mater_param.K_yield=432.2239;
mater_param.n_x=1;

mater_param.c3=0;
mater_param.kappa3=0;
mater_param.eta3=1e14;

mater_param.B=1.3001e-10*100^(y(1)-1);
mater_param.m=4.1913*exp(y(2)-1);
mater_param.D=0.0102*exp(y(3)-1);
mater_param.sigma1=513.3093*exp(y(4)-1);
mater_param.A_nuc=9.1890e-4*10^(y(5)-1);

% mater_param.B=0;
% mater_param.m=1;
% mater_param.D=0;
% mater_param.sigma1=1;
% mater_param.A_nuc=0;

[res_creep_rate, res_stress]=Error_CF_Fortran(mater_param);
res=res_creep_rate+res_stress;
load('bestresult.mat','bestvalue', 'bestparameter');
if(res < bestvalue)
    disp(y);
    bestvalue=res;
    bestparameter=y;
    save('bestresult.mat','bestvalue', 'bestparameter');
    stress_amp_resid=res_stress
    creep_rate_resid=res_creep_rate
   %LCF_resid=res2
end
end