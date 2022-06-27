        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 13 13:32:38 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ONESTEPNAFKRIMPLICIT_EXPLICIT__genmod
          INTERFACE 
            SUBROUTINE ONESTEPNAFKRIMPLICIT_EXPLICIT(EPS_CR,EPS_II,SIGMA&
     &,X,OMEGA,EPS,NEPS_CR,NEPS_II,NOMEGA,K,MU,A,N_NORTON,KAPPA,C,B,M,DT&
     &)
              REAL(KIND=8) :: EPS_CR(3,3)
              REAL(KIND=8) :: EPS_II(3,3)
              REAL(KIND=8) :: SIGMA(3,3)
              REAL(KIND=8) :: X(3,3)
              REAL(KIND=8) :: OMEGA
              REAL(KIND=8) :: EPS(3,3)
              REAL(KIND=8) :: NEPS_CR(3,3)
              REAL(KIND=8) :: NEPS_II(3,3)
              REAL(KIND=8) :: NOMEGA
              REAL(KIND=8) :: K
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: A
              REAL(KIND=8) :: N_NORTON
              REAL(KIND=8) :: KAPPA
              REAL(KIND=8) :: C
              REAL(KIND=8) :: B
              REAL(KIND=8) :: M
              REAL(KIND=8) :: DT
            END SUBROUTINE ONESTEPNAFKRIMPLICIT_EXPLICIT
          END INTERFACE 
        END MODULE ONESTEPNAFKRIMPLICIT_EXPLICIT__genmod
