MODULE EFTDIR_MOD
CONTAINS
SUBROUTINE EFTDIR(PREEL, KFIELDS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRBT, JPRB

USE TPM_DISTR       ,ONLY : D
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE TPM_DIM         ,ONLY : R
USE CUDA_DEVICE_MOD
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IRLEN,ICLEN
INTEGER(KIND=JPIM) :: IPLAN_R2C
REAL(KIND=JPRBT)   :: ZSCAL

integer :: istat

!     ------------------------------------------------------------------

#ifdef UNDEF
BLOCK
INTEGER :: JF, JJ
WRITE (*, *) __FILE__, ':', __LINE__
!$acc serial present (PREEL)
DO JF = 1, SIZE (PREEL, 2)

  PRINT *, "----", JF

  DO JJ = 1, SIZE (PREEL, 1)
    PRINT *, JJ, PREEL (JJ, JF)
  ENDDO
  
ENDDO
!$acc end serial
ENDBLOCK
#endif


BLOCK
INTEGER :: JF, JJ
WRITE (*, *) __FILE__, ':', __LINE__
!$acc serial
DO JF = 1, SIZE (PREEL, 2)

  PRINT *, "----", JF

  DO JJ = 1, SIZE (PREEL, 1)
    IF (MODULO (JJ, 3) == 0) THEN
      PREEL (JJ, JF) = -1._8
    ELSE
      PREEL (JJ, JF) = +1._8
    ENDIF
    PRINT *, JJ, PREEL (JJ, JF)
  ENDDO
  
ENDDO
!$acc end serial
ENDBLOCK


IRLEN=R%NDLON+R%NNOEXTZG
ICLEN=D%NLENGTF/D%NDGL_FS

CALL CREATE_PLAN_FFT (IPLAN_R2C, -1, KN=IRLEN, KLOT=KFIELDS*D%NDGL_FS, &
                    & KISTRIDE=1, KIDIST=ICLEN, KOSTRIDE=1, KODIST=ICLEN/2)
!$acc host_data use_device(PREEL)
CALL EXECUTE_PLAN_FFTC (IPLAN_R2C, -1, PREEL (1, 1))
!$acc end host_data

#ifdef UNDEF
BLOCK
INTEGER :: JF, JJ
WRITE (*, *) __FILE__, ':', __LINE__
!$acc serial present (PREEL)
DO JF = 1, SIZE (PREEL, 2)

  PRINT *, "----", JF

  DO JJ = 1, SIZE (PREEL, 1)
    PRINT *, JJ, PREEL (JJ, JF)
  ENDDO
  
ENDDO
!$acc end serial
ENDBLOCK
#endif

istat = cuda_Synchronize()

ZSCAL = 1._JPRB / REAL (R%NDLON, JPRB)

!$acc kernels present (PREEL) copyin (ZSCAL)
PREEL = ZSCAL * PREEL
!$acc end kernels

BLOCK
INTEGER :: JF, JJ
WRITE (*, *) __FILE__, ':', __LINE__
!$acc serial present (PREEL)
DO JF = 1, SIZE (PREEL, 2)

  PRINT *, "----", JF

  DO JJ = 1, SIZE (PREEL, 1)
    PRINT *, JJ, PREEL (JJ, JF)
  ENDDO
  
ENDDO
!$acc end serial
ENDBLOCK

!     ------------------------------------------------------------------

END SUBROUTINE EFTDIR
END MODULE EFTDIR_MOD
