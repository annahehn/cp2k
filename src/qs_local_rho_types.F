!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE qs_local_rho_types

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE memory_utilities,                ONLY: reallocate
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_rho0_types,                   ONLY: deallocate_rho0_atom,&
                                              deallocate_rho0_mpole,&
                                              rho0_atom_type,&
                                              rho0_mpole_type
   USE qs_rho_atom_types,               ONLY: deallocate_rho_atom_set,&
                                              rho_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_local_rho_types'

! *** Define rhoz and local_rho types ***

! **************************************************************************************************
   TYPE rhoz_type
      REAL(dp)                             ::  one_atom = -1.0_dp
      REAL(dp), DIMENSION(:), POINTER      ::  r_coef => NULL()
      REAL(dp), DIMENSION(:), POINTER      ::  dr_coef => NULL()
      REAL(dp), DIMENSION(:), POINTER      ::  vr_coef => NULL()
   END TYPE rhoz_type

! **************************************************************************************************
   TYPE local_rho_type
      TYPE(rho_atom_type), DIMENSION(:), POINTER            :: rho_atom_set => NULL()
      TYPE(rho0_mpole_type), POINTER                        :: rho0_mpole => NULL()
      TYPE(rho0_atom_type), DIMENSION(:), POINTER           :: rho0_atom_set => NULL()
      TYPE(rhoz_type), DIMENSION(:), POINTER                :: rhoz_set => NULL()
      REAL(dp)                                              :: rhoz_tot = -1.0_dp
   END TYPE local_rho_type

! Public Types
   PUBLIC :: local_rho_type, rhoz_type

! Public Subroutine
   PUBLIC :: allocate_rhoz, calculate_rhoz, &
             get_local_rho, local_rho_set_create, &
             local_rho_set_release, set_local_rho

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rhoz_set ...
!> \param nkind ...
! **************************************************************************************************
   SUBROUTINE allocate_rhoz(rhoz_set, nkind)

      TYPE(rhoz_type), DIMENSION(:), POINTER             :: rhoz_set
      INTEGER                                            :: nkind

      INTEGER                                            :: ikind

      IF (ASSOCIATED(rhoz_set)) THEN
         CALL deallocate_rhoz(rhoz_set)
      END IF

      ALLOCATE (rhoz_set(nkind))

      DO ikind = 1, nkind
         NULLIFY (rhoz_set(ikind)%r_coef)
         NULLIFY (rhoz_set(ikind)%dr_coef)
         NULLIFY (rhoz_set(ikind)%vr_coef)
      END DO

   END SUBROUTINE allocate_rhoz

! **************************************************************************************************
!> \brief ...
!> \param rhoz ...
!> \param grid_atom ...
!> \param alpha ...
!> \param zeff ...
!> \param natom ...
!> \param rhoz_tot ...
!> \param harmonics ...
! **************************************************************************************************
   SUBROUTINE calculate_rhoz(rhoz, grid_atom, alpha, zeff, natom, rhoz_tot, harmonics)

      TYPE(rhoz_type)                                    :: rhoz
      TYPE(grid_atom_type)                               :: grid_atom
      REAL(dp), INTENT(IN)                               :: alpha
      REAL(dp)                                           :: zeff
      INTEGER                                            :: natom
      REAL(dp), INTENT(INOUT)                            :: rhoz_tot
      TYPE(harmonics_atom_type)                          :: harmonics

      INTEGER                                            :: ir, na, nr
      REAL(dp)                                           :: c1, c2, c3, prefactor1, prefactor2, &
                                                            prefactor3, sum

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      CALL reallocate(rhoz%r_coef, 1, nr)
      CALL reallocate(rhoz%dr_coef, 1, nr)
      CALL reallocate(rhoz%vr_coef, 1, nr)

      c1 = alpha/pi
      c2 = c1*c1*c1*fourpi
      c3 = SQRT(alpha)
      prefactor1 = zeff*SQRT(c2)
      prefactor2 = -2.0_dp*alpha
      prefactor3 = -zeff*SQRT(fourpi)

      sum = 0.0_dp
      DO ir = 1, nr
         c1 = -alpha*grid_atom%rad2(ir)
         rhoz%r_coef(ir) = -EXP(c1)*prefactor1
         IF (ABS(rhoz%r_coef(ir)) < 1.0E-30_dp) THEN
            rhoz%r_coef(ir) = 0.0_dp
            rhoz%dr_coef(ir) = 0.0_dp
         ELSE
            rhoz%dr_coef(ir) = prefactor2*rhoz%r_coef(ir)
         END IF
         rhoz%vr_coef(ir) = prefactor3*erf(grid_atom%rad(ir)*c3)/grid_atom%rad(ir)
         sum = sum + rhoz%r_coef(ir)*grid_atom%wr(ir)
      END DO
      rhoz%one_atom = sum*harmonics%slm_int(1)
      rhoz_tot = rhoz_tot + natom*rhoz%one_atom

   END SUBROUTINE calculate_rhoz

! **************************************************************************************************
!> \brief ...
!> \param rhoz_set ...
! **************************************************************************************************
   SUBROUTINE deallocate_rhoz(rhoz_set)

      TYPE(rhoz_type), DIMENSION(:), POINTER             :: rhoz_set

      INTEGER                                            :: ikind, nkind

      nkind = SIZE(rhoz_set)

      DO ikind = 1, nkind
         DEALLOCATE (rhoz_set(ikind)%r_coef)
         DEALLOCATE (rhoz_set(ikind)%dr_coef)
         DEALLOCATE (rhoz_set(ikind)%vr_coef)
      END DO

      DEALLOCATE (rhoz_set)

   END SUBROUTINE deallocate_rhoz

! **************************************************************************************************
!> \brief ...
!> \param local_rho_set ...
!> \param rho_atom_set ...
!> \param rho0_atom_set ...
!> \param rho0_mpole ...
!> \param rhoz_set ...
! **************************************************************************************************
   SUBROUTINE get_local_rho(local_rho_set, rho_atom_set, rho0_atom_set, rho0_mpole, rhoz_set)

      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(rho_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_atom_set
      TYPE(rho0_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho0_atom_set
      TYPE(rho0_mpole_type), OPTIONAL, POINTER           :: rho0_mpole
      TYPE(rhoz_type), DIMENSION(:), OPTIONAL, POINTER   :: rhoz_set

      IF (PRESENT(rho_atom_set)) rho_atom_set => local_rho_set%rho_atom_set
      IF (PRESENT(rho0_atom_set)) rho0_atom_set => local_rho_set%rho0_atom_set
      IF (PRESENT(rho0_mpole)) rho0_mpole => local_rho_set%rho0_mpole
      IF (PRESENT(rhoz_set)) rhoz_set => local_rho_set%rhoz_set

   END SUBROUTINE get_local_rho

! **************************************************************************************************
!> \brief ...
!> \param local_rho_set ...
! **************************************************************************************************
   SUBROUTINE local_rho_set_create(local_rho_set)

      TYPE(local_rho_type), POINTER                      :: local_rho_set

      ALLOCATE (local_rho_set)

      NULLIFY (local_rho_set%rho_atom_set)
      NULLIFY (local_rho_set%rho0_atom_set)
      NULLIFY (local_rho_set%rho0_mpole)
      NULLIFY (local_rho_set%rhoz_set)

      local_rho_set%rhoz_tot = 0.0_dp

   END SUBROUTINE local_rho_set_create

! **************************************************************************************************
!> \brief ...
!> \param local_rho_set ...
! **************************************************************************************************
   SUBROUTINE local_rho_set_release(local_rho_set)

      TYPE(local_rho_type), POINTER                      :: local_rho_set

      IF (ASSOCIATED(local_rho_set)) THEN
         IF (ASSOCIATED(local_rho_set%rho_atom_set)) THEN
            CALL deallocate_rho_atom_set(local_rho_set%rho_atom_set)
         END IF

         IF (ASSOCIATED(local_rho_set%rho0_atom_set)) THEN
            CALL deallocate_rho0_atom(local_rho_set%rho0_atom_set)
         END IF

         IF (ASSOCIATED(local_rho_set%rho0_mpole)) THEN
            CALL deallocate_rho0_mpole(local_rho_set%rho0_mpole)
         END IF

         IF (ASSOCIATED(local_rho_set%rhoz_set)) THEN
            CALL deallocate_rhoz(local_rho_set%rhoz_set)
         END IF

         DEALLOCATE (local_rho_set)
      END IF

   END SUBROUTINE local_rho_set_release

! **************************************************************************************************
!> \brief ...
!> \param local_rho_set ...
!> \param rho_atom_set ...
!> \param rho0_atom_set ...
!> \param rho0_mpole ...
!> \param rhoz_set ...
! **************************************************************************************************
   SUBROUTINE set_local_rho(local_rho_set, rho_atom_set, rho0_atom_set, rho0_mpole, &
                            rhoz_set)

      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(rho_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_atom_set
      TYPE(rho0_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho0_atom_set
      TYPE(rho0_mpole_type), OPTIONAL, POINTER           :: rho0_mpole
      TYPE(rhoz_type), DIMENSION(:), OPTIONAL, POINTER   :: rhoz_set

      IF (PRESENT(rho_atom_set)) THEN
         IF (ASSOCIATED(local_rho_set%rho_atom_set)) THEN
            CALL deallocate_rho_atom_set(local_rho_set%rho_atom_set)
         END IF
         local_rho_set%rho_atom_set => rho_atom_set
      END IF

      IF (PRESENT(rho0_atom_set)) THEN
         IF (ASSOCIATED(local_rho_set%rho0_atom_set)) THEN
            CALL deallocate_rho0_atom(local_rho_set%rho0_atom_set)
         END IF
         local_rho_set%rho0_atom_set => rho0_atom_set
      END IF

      IF (PRESENT(rho0_mpole)) THEN
         IF (ASSOCIATED(local_rho_set%rho0_mpole)) THEN
            CALL deallocate_rho0_mpole(local_rho_set%rho0_mpole)
         END IF
         local_rho_set%rho0_mpole => rho0_mpole
      END IF

      IF (PRESENT(rhoz_set)) THEN
         IF (ASSOCIATED(local_rho_set%rhoz_set)) THEN
            CALL deallocate_rhoz(local_rho_set%rhoz_set)
         END IF
         local_rho_set%rhoz_set => rhoz_set
      END IF

   END SUBROUTINE set_local_rho

END MODULE qs_local_rho_types

