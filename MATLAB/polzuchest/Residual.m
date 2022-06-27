% Функция невязки, которая по параметрам даёт вектор невязки
function [residu]=Residual(x)
load('C:\Users\Viktor\Desktop\ED63.txt');
sizeexp=size(ED63);
Nexp63=sizeexp(1);
e_exp_63=ED63(:,1);
sigma_exp_63=ED63(:,2);
load('C:\Users\Viktor\Desktop\ED53.txt');
sizeexp=size(ED53);
Nexp53=sizeexp(1);
e_exp_53=ED53(:,1);
sigma_exp_53=ED53(:,2);
load('C:\Users\Viktor\Desktop\ED44.txt');
sizeexp=size(ED44);
Nexp44=sizeexp(1);
e_exp_44=ED44(:,1);
sigma_exp_44=ED44(:,2);

% Перевод из безразмерных параметров x в физические параметры, такие как gamma, beta, c1, c2, c3
N=300;
gamma=40000*(-0.22)*x(1);
beta=max(0.1, 100*0.31*x(2));
c(1)=max(0, 120000*1.95*x(3));
c(2)=max(0, 40000*3.52*x(4));
c(3)=20000*0;
kappa(1)=max(0, x(5)*1.35*1/50);
kappa(2)=max(0, x(6)*0.39*1/400);
kappa(3)=1;
K=max(0, x(7)*440);
[e_num_63, sigma_num_63] = SimulateCyclicTension(gamma, beta, c, kappa, K, 0.0063);
[e_num_53, sigma_num_53] = SimulateCyclicTension(gamma, beta, c, kappa, K, 0.0053);
[e_num_44, sigma_num_44] = SimulateCyclicTension(gamma, beta, c, kappa, K, 0.0044);
plot(e_num_63(201:300), sigma_num_63(201:300), 'k');
hold on;
plot(e_num_53(201:300), sigma_num_53(201:300), 'g');
hold on;
plot(e_num_44(201:300), sigma_num_44(201:300), 'b');
hold on;
plot(e_exp_63, sigma_exp_63, 'r');
plot(e_exp_53, sigma_exp_53, 'y');
plot(e_exp_44, sigma_exp_44, 'm');
hold off;
pause(0.1);
eexpacc_63(1)=e_exp_63(1);
for j=2:Nexp63
    eexpacc_63(j)=eexpacc_63(j-1) + abs(e_exp_63(j)-e_exp_63(j-1));
end
enumacc_63(201)=e_num_63(201);
for j=202:N
    enumacc_63(j)=enumacc_63(j-1) + abs(e_num_63(j)-e_num_63(j-1));
end
for j=1:Nexp63
    n=201;
    mn=n;
    r=abs(eexpacc_63(j)-enumacc_63(n));
    for n=202:300
       if r>abs(eexpacc_63(j)-enumacc_63(n))
            r= abs(eexpacc_63(j)-enumacc_63(n));
            mn=n;
       end
    end
    residu(j)=sigma_exp_63(j)-sigma_num_63(mn);
end
eexpacc_53(1)=e_exp_53(1);
for j=2:Nexp53
    eexpacc_53(j)=eexpacc_53(j-1) + abs(e_exp_53(j)-e_exp_53(j-1));
end
enumacc_53(201)=e_num_53(201);
for j=202:N
    enumacc_53(j)=enumacc_53(j-1) + abs(e_num_53(j)-e_num_53(j-1));
end
for j=1:Nexp53
    n=201;
    mn=n;
    r=abs(eexpacc_53(j)-enumacc_53(n));
    for n=202:300
       if r>abs(eexpacc_53(j)-enumacc_53(n))
            r= abs(eexpacc_53(j)-enumacc_53(n));
            mn=n;
       end
    end
    residu(j+Nexp63)=sigma_exp_53(j)-sigma_num_53(mn);
end
eexpacc_44(1)=e_exp_44(1);
for j=2:Nexp44
    eexpacc_44(j)=eexpacc_44(j-1) + abs(e_exp_44(j)-e_exp_44(j-1));
end
enumacc_44(201)=e_num_44(201);
for j=202:N
    enumacc_44(j)=enumacc_44(j-1) + abs(e_num_44(j)-e_num_44(j-1));
end
for j=1:Nexp44
    n=201;
    mn=n;
    r=abs(eexpacc_44(j)-enumacc_44(n));
    for n=202:300
       if r>abs(eexpacc_44(j)-enumacc_44(n))
            r= abs(eexpacc_44(j)-enumacc_44(n));
            mn=n;
       end
    end
    residu(j+Nexp63+Nexp53)=sigma_exp_44(j)-sigma_num_44(mn);
end