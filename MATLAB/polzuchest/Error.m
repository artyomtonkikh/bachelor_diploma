function [res] = Error(x)
A=x(1)*8.4489e-15;
n_Norton=x(2)*2.9054;
c=x(3)*15633;
kappa=x(4)*0.01265;
k=80820.63;
mu=58525.28;
res=Error244(A,n_Norton, c, kappa,k,mu)+1.5*Error266(A,n_Norton, c, kappa,k,mu);
res=res+0;
end

