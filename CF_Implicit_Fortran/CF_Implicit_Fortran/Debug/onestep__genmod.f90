        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 19 19:39:51 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ONESTEP__genmod
          INTERFACE 
            SUBROUTINE ONESTEP(EPS_I,EPS_II1,EPS_II2,SIGMA,F,OMEGA,EPS, &
     &NEPS_I,NEPS_II1,NEPS_II2,NOMEGA,NF,K_INIT,MU_INIT,C1_INIT,C2_INIT,&
     &KAPPA1_INIT,KAPPA2_INIT,ETA1,ETA2,N_NORTON,A_NORTON_INIT,N_X,     &
     &K_YIELD_INIT,D,SIGMA1,ETA,A_NUC,B,M,DT)
              REAL(KIND=8) :: EPS_I(3,3)
              REAL(KIND=8) :: EPS_II1(3,3)
              REAL(KIND=8) :: EPS_II2(3,3)
              REAL(KIND=8) :: SIGMA(3,3)
              REAL(KIND=8) :: F
              REAL(KIND=8) :: OMEGA
              REAL(KIND=8) :: EPS(3,3)
              REAL(KIND=8) :: NEPS_I(3,3)
              REAL(KIND=8) :: NEPS_II1(3,3)
              REAL(KIND=8) :: NEPS_II2(3,3)
              REAL(KIND=8) :: NOMEGA
              REAL(KIND=8) :: NF
              REAL(KIND=8) :: K_INIT
              REAL(KIND=8) :: MU_INIT
              REAL(KIND=8) :: C1_INIT
              REAL(KIND=8) :: C2_INIT
              REAL(KIND=8) :: KAPPA1_INIT
              REAL(KIND=8) :: KAPPA2_INIT
              REAL(KIND=8) :: ETA1
              REAL(KIND=8) :: ETA2
              REAL(KIND=8) :: N_NORTON
              REAL(KIND=8) :: A_NORTON_INIT
              REAL(KIND=8) :: N_X
              REAL(KIND=8) :: K_YIELD_INIT
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
