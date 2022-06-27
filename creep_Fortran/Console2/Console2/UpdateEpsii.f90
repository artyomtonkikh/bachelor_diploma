subroutine UpdateEpsii(eps_ii, neps_ii, eps_cr, ksi, kappa, c)
    !input: neps_ii, eps_cr, ksi, kappa, c
    !output: eps_ii
    real*8 :: eps_ii(3,3), neps_ii(3,3), eps_cr(3,3)
    real*8 :: ksi, kappa, c
    eps_ii=(neps_ii+ksi*kappa*c*eps_cr)/(1.0d0+ksi*kappa*c)
end subroutine UpdateEpsii