MODULE ELEDIR_MOD
CONTAINS
SUBROUTINE ELEDIR(KFC,KLED2,PFFA)

!**** *ELEDIR* - Direct meridional transform.

!     Purpose.
!     --------
!        Direct meridional tranform of state variables.

!**   Interface.
!     ----------
!        CALL ELEDIR(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSIA - symmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              POA1 -  spectral
!                              fields for zonal wavenumber KM
!                              PLEPO - Legendre polonomials

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Reference.
!     ----------

!     Author.
!     -------

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DISTR       ,ONLY : D, D_NUMP
USE TPM_DIM         ,ONLY : R
USE TPMALD_FFT      ,ONLY : TALD
USE TPMALD_DIM      ,ONLY : RALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFC,KLED2
REAL(KIND=JPRB) ,   INTENT(INOUT)  :: PFFA(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: KMLOC, JF, JJ
REAL (KIND=JPRB)   :: ZSCAL
!     ------------------------------------------------------------------

!*       1.       PERFORM FOURIER TRANFORM.
!                 --------------------------

IRLEN=R%NDGL+R%NNOEXTZG
ICLEN=RALD%NDGLSUR+R%NNOEXTZG

CALL CREATE_PLAN_FFT (IPLAN_R2C, -1, KN=IRLEN, KLOT=SIZE (PFFA,2)*SIZE (PFFA, 3), &
                    & KISTRIDE=1, KIDIST=ICLEN, KOSTRIDE=1, KODIST=ICLEN/2)

!$acc host_data use_device (PFFA) 
CALL EXECUTE_PLAN_FFTC(IPLAN_R2C, -1, PFFA (1, 1, 1))
!$acc end host_data

ZSCAL = 1._JPRB / REAL (IRLEN, JPRB)

!$acc parallel loop collapse (3) copyin (D_NUMP, KFC, ICLEN, ZSCAL) present (PFFA)
DO KMLOC = 1, D_NUMP
  DO JF = 1, KFC
    DO JJ = 1, ICLEN
      PFFA (JJ, JF, KMLOC) = PFFA (JJ, JF, KMLOC) * ZSCAL
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

!     ------------------------------------------------------------------

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD
