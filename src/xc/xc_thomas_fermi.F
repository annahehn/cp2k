!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the Thomas-Fermi kinetic energy functional
!> \note
!>      Order of derivatives is: LDA 0; 1; 2; 3;
!>                               LSD 0; a  b; aa bb; aaa bbb;
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \author JGH (18.02.2002)
! **************************************************************************************************
MODULE xc_thomas_fermi
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE kinds,                           ONLY: dp
   USE xc_derivative_desc,              ONLY: deriv_rho,&
                                              deriv_rhoa,&
                                              deriv_rhob
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_get_derivative
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_functionals_utilities,        ONLY: set_util
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
   REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                               f23 = 2.0_dp*f13, &
                               f43 = 4.0_dp*f13, &
                               f53 = 5.0_dp*f13

   PUBLIC :: thomas_fermi_info, thomas_fermi_lda_eval, thomas_fermi_lsd_eval

   REAL(KIND=dp) :: cf, flda, flsd
   REAL(KIND=dp) :: eps_rho

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_thomas_fermi'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cutoff ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_init(cutoff)

      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      eps_rho = cutoff
      CALL set_util(cutoff)

      cf = 0.3_dp*(3.0_dp*pi*pi)**f23
      flda = cf
      flsd = flda*2.0_dp**f23

   END SUBROUTINE thomas_fermi_init

! **************************************************************************************************
!> \brief ...
!> \param lsd ...
!> \param reference ...
!> \param shortform ...
!> \param needs ...
!> \param max_deriv ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_info(lsd, reference, shortform, needs, max_deriv)
      LOGICAL, INTENT(in)                                :: lsd
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "Thomas-Fermi kinetic energy functional: see Parr and Yang"
         IF (.NOT. lsd) THEN
            IF (LEN_TRIM(reference) + 6 < LEN(reference)) THEN
               reference(LEN_TRIM(reference):LEN_TRIM(reference) + 6) = ' {LDA}'
            END IF
         END IF
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "Thomas-Fermi kinetic energy functional"
         IF (.NOT. lsd) THEN
            IF (LEN_TRIM(shortform) + 6 < LEN(shortform)) THEN
               shortform(LEN_TRIM(shortform):LEN_TRIM(shortform) + 6) = ' {LDA}'
            END IF
         END IF
      END IF
      IF (PRESENT(needs)) THEN
         IF (lsd) THEN
            needs%rho_spin = .TRUE.
            needs%rho_spin_1_3 = .TRUE.
         ELSE
            needs%rho = .TRUE.
            needs%rho_1_3 = .TRUE.
         END IF
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 3

   END SUBROUTINE thomas_fermi_info

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lda_eval(rho_set, deriv_set, order)
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: order

      CHARACTER(len=*), PARAMETER :: routineN = 'thomas_fermi_lda_eval'

      INTEGER                                            :: handle, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: epsilon_rho
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: e_0, e_rho, e_rho_rho, e_rho_rho_rho, &
                                                            r13, rho
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)

      CALL xc_rho_set_get(rho_set, rho_1_3=r13, rho=rho, &
                          local_bounds=bo, rho_cutoff=epsilon_rho)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)
      CALL thomas_fermi_init(epsilon_rho)

      IF (order >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)

         CALL thomas_fermi_lda_0(rho, r13, e_0, npoints)
      END IF
      IF (order >= 1 .OR. order == -1) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho)

         CALL thomas_fermi_lda_1(rho, r13, e_rho, npoints)
      END IF
      IF (order >= 2 .OR. order == -2) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)

         CALL thomas_fermi_lda_2(rho, r13, e_rho_rho, npoints)
      END IF
      IF (order >= 3 .OR. order == -3) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)

         CALL thomas_fermi_lda_3(rho, r13, e_rho_rho_rho, npoints)
      END IF
      IF (order > 3 .OR. order < -3) THEN
         CPABORT("derivatives bigger than 3 not implemented")
      END IF
      CALL timestop(handle)
   END SUBROUTINE thomas_fermi_lda_eval

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lsd_eval(rho_set, deriv_set, order)
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: order

      CHARACTER(len=*), PARAMETER :: routineN = 'thomas_fermi_lsd_eval'
      INTEGER, DIMENSION(2), PARAMETER :: rho_spin_name = [deriv_rhoa, deriv_rhob]

      INTEGER                                            :: handle, i, ispin, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: epsilon_rho
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: e_0, e_rho, e_rho_rho, e_rho_rho_rho
      TYPE(cp_3d_r_cp_type), DIMENSION(2)                :: rho, rho_1_3
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)
      NULLIFY (deriv)
      DO i = 1, 2
         NULLIFY (rho(i)%array, rho_1_3(i)%array)
      END DO

      CALL xc_rho_set_get(rho_set, rhoa_1_3=rho_1_3(1)%array, &
                          rhob_1_3=rho_1_3(2)%array, rhoa=rho(1)%array, &
                          rhob=rho(2)%array, &
                          rho_cutoff=epsilon_rho, &
                          local_bounds=bo)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)
      CALL thomas_fermi_init(epsilon_rho)

      DO ispin = 1, 2
         IF (order >= 0) THEN
            deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_0)

            CALL thomas_fermi_lsd_0(rho(ispin)%array, rho_1_3(ispin)%array, &
                                    e_0, npoints)
         END IF
         IF (order >= 1 .OR. order == -1) THEN
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho)

            CALL thomas_fermi_lsd_1(rho(ispin)%array, rho_1_3(ispin)%array, &
                                    e_rho, npoints)
         END IF
         IF (order >= 2 .OR. order == -2) THEN
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        rho_spin_name(ispin)], allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)

            CALL thomas_fermi_lsd_2(rho(ispin)%array, rho_1_3(ispin)%array, &
                                    e_rho_rho, npoints)
         END IF
         IF (order >= 3 .OR. order == -3) THEN
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        rho_spin_name(ispin), rho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)

            CALL thomas_fermi_lsd_3(rho(ispin)%array, rho_1_3(ispin)%array, &
                                    e_rho_rho_rho, npoints)
         END IF
         IF (order > 3 .OR. order < -3) THEN
            CPABORT("derivatives bigger than 3 not implemented")
         END IF
      END DO
      CALL timestop(handle)
   END SUBROUTINE thomas_fermi_lsd_eval

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param e_0 ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lda_0(rho, r13, e_0, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, r13
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_0,flda,r13)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_0(ip) = e_0(ip) + flda*r13(ip)*r13(ip)*rho(ip)

         END IF

      END DO

   END SUBROUTINE thomas_fermi_lda_0

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param e_rho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lda_1(rho, r13, e_rho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, r13
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f53*flda

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE) &
!$OMP SHARED(npoints,rho,eps_rho,e_rho,f,r13)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_rho(ip) = e_rho(ip) + f*r13(ip)*r13(ip)

         END IF

      END DO

   END SUBROUTINE thomas_fermi_lda_1

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param e_rho_rho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lda_2(rho, r13, e_rho_rho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, r13
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f23*f53*flda

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_rho_rho,f,r13)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_rho_rho(ip) = e_rho_rho(ip) + f/r13(ip)

         END IF

      END DO

   END SUBROUTINE thomas_fermi_lda_2

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param e_rho_rho_rho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lda_3(rho, r13, e_rho_rho_rho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, r13
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho_rho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = -f13*f23*f53*flda

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_rho_rho_rho,f,r13)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_rho_rho_rho(ip) = e_rho_rho_rho(ip) + f/(r13(ip)*rho(ip))

         END IF

      END DO

   END SUBROUTINE thomas_fermi_lda_3

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param r13a ...
!> \param e_0 ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lsd_0(rhoa, r13a, e_0, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, r13a
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,e_0,flsd,r13a)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_0(ip) = e_0(ip) + flsd*r13a(ip)*r13a(ip)*rhoa(ip)
         END IF

      END DO

   END SUBROUTINE thomas_fermi_lsd_0

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param r13a ...
!> \param e_rho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lsd_1(rhoa, r13a, e_rho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, r13a
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f53*flsd

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,e_rho,f,r13a)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_rho(ip) = e_rho(ip) + f*r13a(ip)*r13a(ip)
         END IF

      END DO

   END SUBROUTINE thomas_fermi_lsd_1

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param r13a ...
!> \param e_rho_rho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lsd_2(rhoa, r13a, e_rho_rho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, r13a
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f23*f53*flsd

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,e_rho_rho,f,r13a)

      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_rho_rho(ip) = e_rho_rho(ip) + f/r13a(ip)
         END IF

      END DO

   END SUBROUTINE thomas_fermi_lsd_2

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param r13a ...
!> \param e_rho_rho_rho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE thomas_fermi_lsd_3(rhoa, r13a, e_rho_rho_rho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, r13a
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho_rho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = -f13*f23*f53*flsd

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,e_rho_rho_rho,f,r13a)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_rho_rho_rho(ip) = e_rho_rho_rho(ip) + f/(r13a(ip)*rhoa(ip))
         END IF

      END DO

   END SUBROUTINE thomas_fermi_lsd_3

END MODULE xc_thomas_fermi

