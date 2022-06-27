n=0;
if(n==0)
    sigma_eps_data=load('C:\Users\artyo\Documents\Visual Studio 2015\Projects\CF_Implicit_Fortran\CF_Implicit_Fortran\sigma_eps.txt');
else
    sigma_eps_data=load('C:\Users\artyo\Documents\Visual Studio 2015\Projects\CF_Explicit_Fortran\CF_Explicit_Fortran\sigma_eps.txt');
end
creep_rate_data=load('C:\Users\artyo\Documents\Visual Studio 2015\Projects\CF_Implicit_Fortran\CF_Implicit_Fortran\creep_rate.txt');
sigma_amp_data=load('C:\Users\artyo\Documents\Visual Studio 2015\Projects\CF_Implicit_Fortran\CF_Implicit_Fortran\stress_amplitude.txt');
cycles=creep_rate_data(:,1);
creep_rate=creep_rate_data(:,2);
stress_amplitude=sigma_amp_data(:,2);
eps=sigma_eps_data(:,1);
sigma=sigma_eps_data(:,2);
time=sigma_eps_data(:,3);
figure;
plot(eps,sigma);
xlabel('eps');
ylabel('stress, MPa');
grid on;
figure;
plot(time,sigma);
xlabel('t');
ylabel('stress, MPa');
grid on;
figure;
plot(time,eps);
xlabel('t');
ylabel('eps, MPa');
grid on;
figure;
semilogx(cycles,stress_amplitude);
xlabel('N cycles');
ylabel('stress amplitude, MPa');
grid on;
figure;
loglog(cycles,creep_rate);
xlabel('N cycles');
ylabel('creep rate, 1/h');
grid on;