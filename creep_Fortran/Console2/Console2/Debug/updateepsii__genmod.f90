        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 06 15:46:03 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATEEPSII__genmod
          INTERFACE 
            SUBROUTINE UPDATEEPSII(EPS_II,NEPS_II,EPS_CR,KSI,KAPPA,C)
              REAL(KIND=8) :: EPS_II(3,3)
              REAL(KIND=8) :: NEPS_II(3,3)
              REAL(KIND=8) :: EPS_CR(3,3)
              REAL(KIND=8) :: KSI
              REAL(KIND=8) :: KAPPA
              REAL(KIND=8) :: C
            END SUBROUTINE UPDATEEPSII
          END INTERFACE 
        END MODULE UPDATEEPSII__genmod
