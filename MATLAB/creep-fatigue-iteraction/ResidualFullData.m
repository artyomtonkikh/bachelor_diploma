function [resid]=ResidualFullData(y)
disp(y);
%идеальные первые и вторые графики
% mater_param.A_Norton=7.28e-07*100^(y(1)-1);
% mater_param.n_Norton=3.0376*y(2);
% mater_param.eta1=1.9812e+03*100^(y(3)-1);
% mater_param.eta2=1.2056e+03*100^(y(4)-1);
% mater_param.c1=1.0595e+04*10^(y(5)-1);
% mater_param.c2=3.5388e+04*10^(y(6)-1);
% mater_param.kappa1=0.0083*exp(y(7)-1);
% mater_param.kappa2=0.0021*exp(y(8)-1);
% mater_param.K_yield=382.9207*10^(y(9)-1);
% mater_param.n_x=1+110*abs(y(10)-1);

%тоже хорошие первые и вторые графики
% mater_param.A_Norton=2.7666e-12*100^(y(1)-1);
% mater_param.n_Norton=4.9177*y(2);
% mater_param.eta1=4.3724e+08*100^(y(3)-1);
% mater_param.eta2=263660000*100^(y(4)-1);
% mater_param.c1=12871*100^(y(5)-1);
% mater_param.c2=12918*100^(y(6)-1);
% mater_param.kappa1=0.0025*exp(y(7)-1);
% mater_param.kappa2=0.0018*exp(y(8)-1);
% mater_param.K_yield=363*abs(y(9));
% mater_param.n_x=3.1954*abs(y(10));

%идеально всё
mater_param.A_Norton=4.9017e-6*100^(y(1)-1);
mater_param.n_Norton=2.5125*y(2);
mater_param.eta1=1.6746e3*100^(y(3)-1);
mater_param.eta2=1.3561e5*100^(y(4)-1);
mater_param.c1=10000*10^(y(5)-1);
mater_param.c2=2000*10^(y(6)-1);
mater_param.kappa1=0.01*exp(y(7)-1);
mater_param.kappa2=0.009*exp(y(8)-1);
mater_param.K_yield=393*10^(y(9)-1);
mater_param.n_x=1;

mater_param.c3=0;
mater_param.kappa3=0;
mater_param.eta3=1e14;
mater_param.B=0;
mater_param.m=1;
mater_param.D=0;
mater_param.sigma1=1;
mater_param.A_nuc=0;
resid1=Error_CF_Fortran_RFD(mater_param);
%resid2=Error_LCF_Fortran_RFD(mater_param);
%resid2=Error_CF_amps_RFD(mater_param);
resid=[resid1];
a=0;
a=a+1;
end