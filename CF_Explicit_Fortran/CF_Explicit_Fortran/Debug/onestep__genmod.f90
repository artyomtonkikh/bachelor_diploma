        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 11 20:49:36 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ONESTEP__genmod
          INTERFACE 
            SUBROUTINE ONESTEP(EPS_CR,EPS_P,EPS_II1,EPS_II2,EPS_II3,    &
     &SIGMA,F,OMEGA,EPS,NEPS_CR,NEPS_P,NEPS_II1,NEPS_II2,NEPS_II3,NOMEGA&
     &,NF,K,MU,C1,C2,C3,KAPPA1,KAPPA2,KAPPA3,ETA1,ETA2,ETA3,N_NORTON,   &
     &A_NORTON,N_X,K_YIELD,D,SIGMA1,ETA,A_NUC,B,M,DT)
              REAL(KIND=8) :: EPS_CR(3,3)
              REAL(KIND=8) :: EPS_P(3,3)
              REAL(KIND=8) :: EPS_II1(3,3)
              REAL(KIND=8) :: EPS_II2(3,3)
              REAL(KIND=8) :: EPS_II3(3,3)
              REAL(KIND=8) :: SIGMA(3,3)
              REAL(KIND=8) :: F
              REAL(KIND=8) :: OMEGA
              REAL(KIND=8) :: EPS(3,3)
              REAL(KIND=8) :: NEPS_CR(3,3)
              REAL(KIND=8) :: NEPS_P(3,3)
              REAL(KIND=8) :: NEPS_II1(3,3)
              REAL(KIND=8) :: NEPS_II2(3,3)
              REAL(KIND=8) :: NEPS_II3(3,3)
              REAL(KIND=8) :: NOMEGA
              REAL(KIND=8) :: NF
              REAL(KIND=8) :: K
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: C1
              REAL(KIND=8) :: C2
              REAL(KIND=8) :: C3
              REAL(KIND=8) :: KAPPA1
              REAL(KIND=8) :: KAPPA2
              REAL(KIND=8) :: KAPPA3
              REAL(KIND=8) :: ETA1
              REAL(KIND=8) :: ETA2
              REAL(KIND=8) :: ETA3
              REAL(KIND=8) :: N_NORTON
              REAL(KIND=8) :: A_NORTON
              REAL(KIND=8) :: N_X
              REAL(KIND=8) :: K_YIELD
              REAL(KIND=8) :: D
              REAL(KIND=8) :: SIGMA1
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: A_NUC
              REAL(KIND=8) :: B
              REAL(KIND=8) :: M
              REAL(KIND=8) :: DT
            END SUBROUTINE ONESTEP
          END INTERFACE 
        END MODULE ONESTEP__genmod
