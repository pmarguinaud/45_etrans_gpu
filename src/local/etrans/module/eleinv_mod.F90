MODULE ELEINV_MOD
CONTAINS
SUBROUTINE ELEINV(KFC,KF_OUT_LT,PFFA)

!**** *LEINV* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINV(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PLEPO - Legendre polonomials for zonal
!                              wavenumber KM (input-c)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINV in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 01-Sep-2015 support for FFTW transforms

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DISTR       ,ONLY : D, D_NUMP
USE TPM_DIM         ,ONLY : R
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FFT      ,ONLY : TALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(INOUT)  :: PFFA(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN
INTEGER(KIND=JPIM) :: KMLOC, JF, JJ
INTEGER(KIND=JPIM) :: IPLAN_C2R
REAL (KIND=JPRB)   :: ZSCAL

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',0,ZHOOK_HANDLE)


IRLEN=R%NDGL+R%NNOEXTZG
ICLEN=RALD%NDGLSUR+R%NNOEXTZG

CALL CREATE_PLAN_FFT (IPLAN_C2R, +1, KN=IRLEN, KLOT=UBOUND (PFFA,2)*D%NUMP, &
                    & KISTRIDE=1, KIDIST=ICLEN/2, KOSTRIDE=1, KODIST=ICLEN)

!$acc host_data use_device (PFFA) 
CALL EXECUTE_PLAN_FFTC(IPLAN_C2R, +1, PFFA (1, 1, 1))
!$acc end host_data

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)

END SUBROUTINE ELEINV
END MODULE ELEINV_MOD
