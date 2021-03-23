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

!     Externals.   MXMAOP - matrix multiply
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - NTMAX instead of NSMAX
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 01-Sep-2015 support for FFTW transforms
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DISTR       ,ONLY : D
USE TPM_DIM         ,ONLY : R
!USE TPM_GEOMETRY
!USE TPM_TRANS
USE TPMALD_FFT      ,ONLY : TALD
#ifdef WITH_FFTW
USE TPM_FFTW     ,ONLY : TW, EXEC_EFFTW
#endif
USE TPMALD_DIM      ,ONLY : RALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFC,KLED2
REAL(KIND=JPRB) ,   INTENT(INOUT)  :: PFFA(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN, IOFF, ITYPE
INTEGER(KIND=JPIM) :: IPLAN_R2C
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time
!     ------------------------------------------------------------------

!*       1.       PERFORM FOURIER TRANFORM.
!                 --------------------------

IRLEN=R%NDGL+R%NNOEXTZG
ICLEN=RALD%NDGLSUR+R%NNOEXTZG

PRINT *, " KFC   = ", KFC
PRINT *, " IRLEN = ", IRLEN
PRINT *, " ICLEN = ", ICLEN
PRINT *, " UBOUND (PFFA) = ", UBOUND (PFFA)

BLOCK
INTEGER :: KMLOC, JF, JJ
WRITE (*, *) __FILE__, ':', __LINE__ 
!$acc serial copyin (D%NUMP,IRLEN) present (PFFA)
PFFA = 0.
DO JJ = 1, IRLEN
  PFFA (JJ, 1, 1) = 1._JPRB
ENDDO
DO JJ = 1, IRLEN
  PFFA (JJ, 2, 1) = 1._JPRB
ENDDO
PFFA (:, :, :) = 999999.
PFFA (:, 1, :) = 1.
PFFA (:, 2, :) = 0.
DO KMLOC = 1, D%NUMP
  DO JJ = 1, ICLEN
    PRINT *, KMLOC, JJ, PFFA (JJ, 1, KMLOC), PFFA (JJ, 2, KMLOC)
  ENDDO
ENDDO
!$acc end serial
ENDBLOCK


!CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,IRLEN,KFC*D%NUMP,1,ICLEN)
CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,IRLEN,UBOUND (PFFA,2)*D%NUMP,1,ICLEN)

!$acc host_data use_device (PFFA)
CALL EXECUTE_PLAN_FFTC(IPLAN_R2C, -1, PFFA)
!$acc end host_data

BLOCK
INTEGER :: KMLOC, JF, JJ
WRITE (*, *) __FILE__, ':', __LINE__ 
!$acc serial copyin (D%NUMP,IRLEN) present (PFFA)
DO KMLOC = 1, D%NUMP
  DO JJ = 1, ICLEN
    PRINT *, KMLOC, JJ, PFFA (JJ, 1, KMLOC), PFFA (JJ, 2, KMLOC)
  ENDDO
ENDDO
!$acc end serial
ENDBLOCK

CALL ABOR1 ('ELEDIR')

#ifdef UNDEF
IF (KFC>0) THEN
  ITYPE=-1
  IRLEN=R%NDGL+R%NNOEXTZG
  ICLEN=RALD%NDGLSUR+R%NNOEXTZG
  IF( TALD%LFFT992 )THEN
    CALL FFT992(PFFT,TALD%TRIGSE,TALD%NFAXE,1,ICLEN,IRLEN,KFC,ITYPE)
#ifdef WITH_FFTW
  ELSEIF( TW%LFFTW )THEN
    IOFF=1
    CALL EXEC_EFFTW(ITYPE,IRLEN,ICLEN,IOFF,KFC,LL_ALL,PFFT)
#endif
  ELSE
    CALL ABORT_TRANS('ELEDIR_MOD:ELEDIR: NO FFT PACKAGE SELECTED')
  ENDIF
ENDIF
#endif

!     ------------------------------------------------------------------

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD
