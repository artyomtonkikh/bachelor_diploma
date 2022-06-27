function [eps_cr] = UpdateEpscr(eps, neps_cr, neps_ii, ksi,kappa, c, mu, F)
denominator=1+2*mu*ksi/F+c*ksi/F-c^2*ksi^2*kappa/(F*(1+ksi*kappa*c));
eps_cr=(neps_cr+2*mu*ksi/F*dev(eps)+ksi*c/(F*(1+ksi*kappa*c))*neps_ii)/denominator;
end