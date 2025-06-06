!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types containing essential information for running implicit (iterative)
!>        Poisson solver
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! **************************************************************************************************
MODULE ps_implicit_types

   USE dct,                             ONLY: dct_type,&
                                              dct_type_release
   USE dielectric_types,                ONLY: dielectric_release,&
                                              dielectric_type
   USE dirichlet_bc_types,              ONLY: dbc_release,&
                                              dirichlet_bc_p_type
   USE kinds,                           ONLY: dp
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_implicit_types'

   INTEGER, PARAMETER, PUBLIC           :: PERIODIC_BC = 0, &
                                           MIXED_BC = 1, &
                                           MIXED_PERIODIC_BC = 2, &
                                           NEUMANN_BC = 3

   TYPE :: ps_implicit_parameters
      INTEGER                           :: boundary_condition = PERIODIC_BC
      LOGICAL                           :: zero_initial_guess = .FALSE.
      INTEGER                           :: max_iter = 0
      REAL(KIND=dp)                     :: tol = 0.0_dp
      REAL(KIND=dp)                     :: omega = 0.0_dp
      INTEGER                           :: neumann_directions = 0
   END TYPE ps_implicit_parameters

   TYPE :: ps_implicit_type
      LOGICAL                           :: do_dbc_cube = .FALSE.
      TYPE(dielectric_type), POINTER    :: dielectric => NULL()
      TYPE(dct_type)                    :: dct_env = dct_type()
      TYPE(pw_r3d_rs_type), POINTER            :: initial_guess => NULL()
      TYPE(pw_r3d_rs_type), POINTER            :: v_eps => NULL()
      TYPE(pw_r3d_rs_type), POINTER            :: cstr_charge => NULL()
      REAL(dp), DIMENSION(:), ALLOCATABLE :: initial_lambda
      REAL(dp)                          :: ehartree = 0.0_dp
      REAL(dp)                          :: electric_enthalpy = 0.0_dp
      INTEGER                           :: times_called = 0
      TYPE(dirichlet_bc_p_type), &
         DIMENSION(:), ALLOCATABLE      :: contacts
      REAL(dp), DIMENSION(:, :), ALLOCATABLE :: QS
      REAL(dp), DIMENSION(:, :), ALLOCATABLE :: Rinv
      REAL(dp), DIMENSION(:, :), ALLOCATABLE :: B
      REAL(dp), DIMENSION(:, :), ALLOCATABLE :: Bt
      REAL(dp), DIMENSION(:), ALLOCATABLE :: v_D
      REAL(dp), DIMENSION(:), ALLOCATABLE :: osc_frac
      REAL(dp), DIMENSION(:), ALLOCATABLE :: frequency
      REAL(dp), DIMENSION(:), ALLOCATABLE :: phase
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx_1dto3d
   END TYPE ps_implicit_type

   PUBLIC ps_implicit_type, ps_implicit_parameters
   PUBLIC ps_implicit_release

CONTAINS

! **************************************************************************************************
!> \brief   Deallocates ps_implicit
!> \param ps_implicit_env the implicit_env to be deallocated
!> \param ps_implicit_params implicit env parameters
!> \param pw_pool pool of plane-wave grid
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE ps_implicit_release(ps_implicit_env, ps_implicit_params, pw_pool)
      TYPE(ps_implicit_type), POINTER                    :: ps_implicit_env
      TYPE(ps_implicit_parameters), INTENT(INOUT)        :: ps_implicit_params
      TYPE(pw_pool_type), INTENT(IN), OPTIONAL, POINTER  :: pw_pool

      CHARACTER(len=*), PARAMETER :: routineN = 'ps_implicit_release'

      INTEGER                                            :: boundary_condition, handle
      LOGICAL                                            :: can_give_back, do_dbc_cube

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(ps_implicit_env)) THEN
         can_give_back = PRESENT(pw_pool)
         IF (can_give_back) can_give_back = ASSOCIATED(pw_pool)
         do_dbc_cube = ps_implicit_env%do_dbc_cube

         IF (can_give_back) THEN
            CALL pw_pool%give_back_pw(ps_implicit_env%initial_guess)
            CALL pw_pool%give_back_pw(ps_implicit_env%v_eps)
            CALL pw_pool%give_back_pw(ps_implicit_env%cstr_charge)
            DEALLOCATE (ps_implicit_env%initial_guess, ps_implicit_env%v_eps, ps_implicit_env%cstr_charge)
            CALL dbc_release(ps_implicit_env%contacts, do_dbc_cube, pw_pool=pw_pool)
         ELSE
            CALL ps_implicit_env%initial_guess%release()
            CALL ps_implicit_env%v_eps%release()
            IF (ASSOCIATED(ps_implicit_env%cstr_charge)) THEN
               CALL ps_implicit_env%cstr_charge%release()
               DEALLOCATE (ps_implicit_env%cstr_charge)
            END IF
            DEALLOCATE (ps_implicit_env%initial_guess, ps_implicit_env%v_eps)
            CALL dbc_release(ps_implicit_env%contacts, do_dbc_cube)
         END IF

         DEALLOCATE (ps_implicit_env%initial_lambda)
         DEALLOCATE (ps_implicit_env%B)
         DEALLOCATE (ps_implicit_env%Bt)
         DEALLOCATE (ps_implicit_env%QS)
         DEALLOCATE (ps_implicit_env%Rinv)
         DEALLOCATE (ps_implicit_env%v_D)
         DEALLOCATE (ps_implicit_env%osc_frac)
         DEALLOCATE (ps_implicit_env%frequency)
         DEALLOCATE (ps_implicit_env%phase)
         DEALLOCATE (ps_implicit_env%idx_1dto3d)

         CALL dielectric_release(ps_implicit_env%dielectric, pw_pool)
         boundary_condition = ps_implicit_params%boundary_condition
         IF (boundary_condition .EQ. MIXED_BC .OR. boundary_condition .EQ. NEUMANN_BC) THEN
            CALL dct_type_release(ps_implicit_env%dct_env)
         END IF

         DEALLOCATE (ps_implicit_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ps_implicit_release

END MODULE ps_implicit_types
