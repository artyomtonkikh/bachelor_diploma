        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 06 15:11:23 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATEEPSCR__genmod
          INTERFACE 
            SUBROUTINE UPDATEEPSCR(EPS_CR,EPS,NEPS_CR,NEPS_II,KSI,KAPPA,&
     &C,MU,F)
              REAL(KIND=8) :: EPS_CR(3,3)
              REAL(KIND=8) :: EPS(3,3)
              REAL(KIND=8) :: NEPS_CR(3,3)
              REAL(KIND=8) :: NEPS_II(3,3)
              REAL(KIND=8) :: KSI
              REAL(KIND=8) :: KAPPA
              REAL(KIND=8) :: C
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: F
            END SUBROUTINE UPDATEEPSCR
          END INTERFACE 
        END MODULE UPDATEEPSCR__genmod
