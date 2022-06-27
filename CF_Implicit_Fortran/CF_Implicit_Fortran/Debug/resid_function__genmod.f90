        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 16 20:12:58 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RESID_FUNCTION__genmod
          INTERFACE 
            FUNCTION RESID_FUNCTION(F_PAR,XI,DEV_EPS,EPS_I,EPS_II1,     &
     &EPS_II2,NEPS_I,NEPS_II1,NEPS_II2,K,MU,A_NORTON,N_NORTON,KAPPA1,   &
     &KAPPA2,C1,C2,ETA1,ETA2,DT)
              REAL(KIND=8) :: F_PAR
              REAL(KIND=8) :: XI
              REAL(KIND=8) :: DEV_EPS(3,3)
              REAL(KIND=8) :: EPS_I(3,3)
              REAL(KIND=8) :: EPS_II1(3,3)
              REAL(KIND=8) :: EPS_II2(3,3)
              REAL(KIND=8) :: NEPS_I(3,3)
              REAL(KIND=8) :: NEPS_II1(3,3)
              REAL(KIND=8) :: NEPS_II2(3,3)
              REAL(KIND=8) :: K
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: A_NORTON
              REAL(KIND=8) :: N_NORTON
              REAL(KIND=8) :: KAPPA1
              REAL(KIND=8) :: KAPPA2
              REAL(KIND=8) :: C1
              REAL(KIND=8) :: C2
              REAL(KIND=8) :: ETA1
              REAL(KIND=8) :: ETA2
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: RESID_FUNCTION
            END FUNCTION RESID_FUNCTION
          END INTERFACE 
        END MODULE RESID_FUNCTION__genmod
