subroutine UpdateEpscr(eps_cr, eps, neps_cr, neps_ii, ksi, kappa, c, mu, F)
    !input: eps, neps_cr, neps_ii, ksi, kappa, c, mu, F
    !output: eps_cr
    real*8 :: eps(3,3), eps_cr(3,3), neps_cr(3,3), neps_ii(3,3)
    real*8 :: dev_eps(3,3)
    real*8 :: eye(3,3) !единичный тензор
    real*8 :: ksi, kappa, c, mu, F
    real*8 :: denominator
    denominator=1.0d0+2.0d0*mu*ksi/F+c*ksi/F-c**2.0d0*ksi**2.0d0*kappa/(F*(1.0d0+ksi*kappa*c))
    eye=0.0d0
    do i=1,3
        eye(i,i)=1.0d0
    end do
    call dev(dev_eps, eps)
    eps_cr=(neps_cr+2.0d0*mu*ksi/F*dev_eps+ksi*c/(F*(1.0d0+ksi*kappa*c))*neps_ii)/denominator
end subroutine UpdateEpscr