        !COMPILER-GENERATED INTERFACE MODULE: Wed Mar 03 21:28:59 2021
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
