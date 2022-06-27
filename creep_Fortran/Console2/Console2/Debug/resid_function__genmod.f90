        !COMPILER-GENERATED INTERFACE MODULE: Fri Oct 09 23:41:08 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RESID_FUNCTION__genmod
          INTERFACE 
            FUNCTION RESID_FUNCTION(KSI,EPS,NEPS_CR,NEPS_II,K,MU,A,     &
     &N_NORTON,KAPPA,C,DT)
              REAL(KIND=8) :: KSI
              REAL(KIND=8) :: EPS(3,3)
              REAL(KIND=8) :: NEPS_CR(3,3)
              REAL(KIND=8) :: NEPS_II(3,3)
              REAL(KIND=8) :: K
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: A
              REAL(KIND=8) :: N_NORTON
              REAL(KIND=8) :: KAPPA
              REAL(KIND=8) :: C
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: RESID_FUNCTION
            END FUNCTION RESID_FUNCTION
          END INTERFACE 
        END MODULE RESID_FUNCTION__genmod
