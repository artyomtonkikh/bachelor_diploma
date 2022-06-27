function [eps_ii] = UpdateEpsii(neps_ii, eps_cr,ksi,kappa,c)
eps_ii=(neps_ii+ksi*kappa*c*eps_cr)/(1+ksi*kappa*c);
end