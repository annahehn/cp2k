!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE xc_atom

   USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_next,&
                                              cp_sll_xc_deriv_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE xc,                              ONLY: divide_by_norm_drho,&
                                              xc_calc_2nd_deriv_analytical
   USE xc_derivative_desc,              ONLY: &
        deriv_norm_drho, deriv_norm_drhoa, deriv_norm_drhob, deriv_rho, deriv_rhoa, deriv_rhob, &
        deriv_tau, deriv_tau_a, deriv_tau_b
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_get_derivative
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_derivatives,                  ONLY: xc_functionals_eval
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_atom'

   PUBLIC :: vxc_of_r_new, xc_rho_set_atom_update, xc_2nd_deriv_of_r, fill_rho_set

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param xc_fun_section ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param deriv_order ...
!> \param needs ...
!> \param w ...
!> \param lsd ...
!> \param na ...
!> \param nr ...
!> \param exc ...
!> \param vxc ...
!> \param vxg ...
!> \param vtau ...
!> \param energy_only ...
!> \param epr_xc ...
!> \param adiabatic_rescale_factor ...
! **************************************************************************************************
   SUBROUTINE vxc_of_r_new(xc_fun_section, rho_set, deriv_set, deriv_order, needs, w, &
                           lsd, na, nr, exc, vxc, vxg, vtau, &
                           energy_only, epr_xc, adiabatic_rescale_factor)

! This routine updates rho_set by giving to it the rho and drho that are needed.
! Since for the local densities rho1_h and rho1_s local grids are used it is not possible
! to call xc_rho_set_update.
! As input of this routine one gets rho and drho on a one dimensional grid.
! The grid is the angular grid corresponding to a given point ir_pnt on the radial grid.
! The derivatives are calculated on this one dimensional grid, the results are stored in
! exc, vxc(1:na,ir_pnt,ispin), vxg(1:na,ir_pnt,ispin), vxg_cross(1:na,ir_pnt,ispin)
! Afterwords the arrays containing the derivatives are put to zero so that the routine
! can safely be called for the next radial point ir_pnt

      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: deriv_order
      TYPE(xc_rho_cflags_type), INTENT(in)               :: needs
      REAL(dp), DIMENSION(:, :), POINTER                 :: w
      LOGICAL, INTENT(IN)                                :: lsd
      INTEGER, INTENT(in)                                :: na, nr
      REAL(dp)                                           :: exc
      REAL(dp), DIMENSION(:, :, :), POINTER              :: vxc
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: vxg
      REAL(dp), DIMENSION(:, :, :), POINTER              :: vtau
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_only, epr_xc
      REAL(dp), INTENT(IN), OPTIONAL                     :: adiabatic_rescale_factor

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vxc_of_r_new'

      INTEGER                                            :: handle, ia, idir, ir, my_deriv_order
      LOGICAL                                            :: gradient_f, my_epr_xc, my_only_energy
      REAL(dp)                                           :: my_adiabatic_rescale_factor
      REAL(dp), DIMENSION(:, :, :), POINTER              :: deriv_data
      REAL(KIND=dp)                                      :: drho_cutoff
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      CALL timeset(routineN, handle)
      my_only_energy = .FALSE.
      IF (PRESENT(energy_only)) my_only_energy = energy_only

      IF (PRESENT(adiabatic_rescale_factor)) THEN
         my_adiabatic_rescale_factor = adiabatic_rescale_factor
      ELSE
         my_adiabatic_rescale_factor = 1.0_dp
      END IF

      ! needed for the epr routines
      my_epr_xc = .FALSE.
      IF (PRESENT(epr_xc)) my_epr_xc = epr_xc
      my_deriv_order = deriv_order
      IF (my_epr_xc) my_deriv_order = 2

      gradient_f = (needs%drho_spin .OR. needs%norm_drho_spin .OR. &
                    needs%drho .OR. needs%norm_drho)

!  Calculate the derivatives
      CALL xc_functionals_eval(xc_fun_section, &
                               lsd=lsd, &
                               rho_set=rho_set, &
                               deriv_set=deriv_set, &
                               deriv_order=my_deriv_order)

      CALL xc_rho_set_get(rho_set, drho_cutoff=drho_cutoff)

      NULLIFY (deriv_data)

      IF (my_epr_xc) THEN
         ! nabla v_xc (using the vxg arrays)
         ! there's no point doing this when lsd = false
         IF (lsd) THEN
            deriv_att => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhoa])
            IF (ASSOCIATED(deriv_att)) THEN
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
               DO ir = 1, nr
                  DO ia = 1, na
                     DO idir = 1, 3
                        vxg(idir, ia, ir, 1) = rho_set%drhoa(idir)%array(ia, ir, 1)* &
                                               deriv_data(ia, ir, 1)
                     END DO !idir
                  END DO !ia
               END DO !ir
               NULLIFY (deriv_data)
            END IF
            deriv_att => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_rhob])
            IF (ASSOCIATED(deriv_att)) THEN
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
               DO ir = 1, nr
                  DO ia = 1, na
                     DO idir = 1, 3
                        vxg(idir, ia, ir, 2) = rho_set%drhob(idir)%array(ia, ir, 1)* &
                                               deriv_data(ia, ir, 1)
                     END DO !idir
                  END DO !ia
               END DO !ir
               NULLIFY (deriv_data)
            END IF
         END IF
         !  EXC energy ! is that needed for epr?
         deriv_att => xc_dset_get_derivative(deriv_set, [INTEGER::])
         exc = 0.0_dp
         IF (ASSOCIATED(deriv_att)) THEN
            CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
            DO ir = 1, nr
               DO ia = 1, na
                  exc = exc + deriv_data(ia, ir, 1)*w(ia, ir)
               END DO
            END DO
            NULLIFY (deriv_data)
         END IF
      ELSE
!  EXC energy
         deriv_att => xc_dset_get_derivative(deriv_set, [INTEGER::])
         exc = 0.0_dp
         IF (ASSOCIATED(deriv_att)) THEN
            CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
            DO ir = 1, nr
               DO ia = 1, na
                  exc = exc + deriv_data(ia, ir, 1)*w(ia, ir)
               END DO
            END DO
            NULLIFY (deriv_data)
         END IF
         ! Calculate the potential only if needed
         IF (.NOT. my_only_energy) THEN
!  Derivative with respect to the density
            IF (lsd) THEN
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_rhoa])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vxc(:, :, 1) = deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_rhob])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vxc(:, :, 2) = deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_rho])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vxc(:, :, 1) = vxc(:, :, 1) + deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  vxc(:, :, 2) = vxc(:, :, 2) + deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
            ELSE
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_rho])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vxc(:, :, 1) = deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
            END IF ! lsd

!  Derivatives with respect to the gradient
            IF (lsd) THEN
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  DO ir = 1, nr
                     DO ia = 1, na
                        DO idir = 1, 3
                           IF (rho_set%norm_drhoa(ia, ir, 1) > drho_cutoff) THEN
                              vxg(idir, ia, ir, 1) = rho_set%drhoa(idir)%array(ia, ir, 1)* &
                                                     deriv_data(ia, ir, 1)*w(ia, ir)/ &
                                                     rho_set%norm_drhoa(ia, ir, 1)*my_adiabatic_rescale_factor
                           ELSE
                              vxg(idir, ia, ir, 1) = 0.0_dp
                           END IF
                        END DO
                     END DO
                  END DO
                  NULLIFY (deriv_data)
               END IF
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  DO ir = 1, nr
                     DO ia = 1, na
                        DO idir = 1, 3
                           IF (rho_set%norm_drhob(ia, ir, 1) > drho_cutoff) THEN
                              vxg(idir, ia, ir, 2) = rho_set%drhob(idir)%array(ia, ir, 1)* &
                                                     deriv_data(ia, ir, 1)*w(ia, ir)/ &
                                                     rho_set%norm_drhob(ia, ir, 1)*my_adiabatic_rescale_factor
                           ELSE
                              vxg(idir, ia, ir, 2) = 0.0_dp
                           END IF
                        END DO
                     END DO
                  END DO
                  NULLIFY (deriv_data)
               END IF
               ! Cross Terms
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_norm_drho])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  DO ir = 1, nr
                     DO ia = 1, na
                        DO idir = 1, 3
                           IF (rho_set%norm_drho(ia, ir, 1) > drho_cutoff) THEN
                              vxg(idir, ia, ir, 1:2) = &
                                 vxg(idir, ia, ir, 1:2) + ( &
                                 rho_set%drhoa(idir)%array(ia, ir, 1) + &
                                 rho_set%drhob(idir)%array(ia, ir, 1))* &
                                 deriv_data(ia, ir, 1)*w(ia, ir)/rho_set%norm_drho(ia, ir, 1)* &
                                 my_adiabatic_rescale_factor
                           END IF
                        END DO
                     END DO
                  END DO
                  NULLIFY (deriv_data)
               END IF
            ELSE
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_norm_drho])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  DO ir = 1, nr
                     DO ia = 1, na
                        IF (rho_set%norm_drho(ia, ir, 1) > drho_cutoff) THEN
                           DO idir = 1, 3
                              vxg(idir, ia, ir, 1) = rho_set%drho(idir)%array(ia, ir, 1)* &
                                                     deriv_data(ia, ir, 1)*w(ia, ir)/ &
                                                     rho_set%norm_drho(ia, ir, 1)*my_adiabatic_rescale_factor
                           END DO
                        ELSE
                           vxg(1:3, ia, ir, 1) = 0.0_dp
                        END IF
                     END DO
                  END DO
                  NULLIFY (deriv_data)
               END IF
            END IF ! lsd
!  Derivative with respect to tau
            IF (lsd) THEN
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_tau_a])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vtau(:, :, 1) = deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_tau_b])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vtau(:, :, 2) = deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_tau])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vtau(:, :, 1) = vtau(:, :, 1) + deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  vtau(:, :, 2) = vtau(:, :, 2) + deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
            ELSE
               deriv_att => xc_dset_get_derivative(deriv_set, [deriv_tau])
               IF (ASSOCIATED(deriv_att)) THEN
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
                  vtau(:, :, 1) = deriv_data(:, :, 1)*w(:, :)*my_adiabatic_rescale_factor
                  NULLIFY (deriv_data)
               END IF
            END IF ! lsd
         END IF ! only_energy
      END IF ! epr_xc

      CALL timestop(handle)

   END SUBROUTINE vxc_of_r_new

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param rho1_set ...
!> \param xc_section ...
!> \param deriv_set ...
!> \param w ...
!> \param vxc ...
!> \param vxg ...
!> \param do_triplet ...
! **************************************************************************************************
   SUBROUTINE xc_2nd_deriv_of_r(rho_set, rho1_set, xc_section, &
                                deriv_set, w, vxc, vxg, do_triplet)

! As input of this routine one gets rho and drho on a one dimensional grid.
! The grid is the angular grid corresponding to a given point ir on the radial grid.
! The derivatives are calculated on this one dimensional grid, the results are stored in
! vxc(1:na,ir,ispin), vxg(1:na,ir,ispin), vxg_cross(1:na,ir,ispin)
! Afterwords the arrays containing the derivatives are put to zero so that the routine
! can safely be called for the next radial point ir

      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set, rho1_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(xc_derivative_set_type), INTENT(INOUT)        :: deriv_set
      REAL(dp), DIMENSION(:, :), POINTER                 :: w
      REAL(dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER  :: vxc
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: vxg
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_triplet

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'xc_2nd_deriv_of_r'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: lsd
      REAL(dp)                                           :: drho_cutoff, my_fac_triplet
      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_pw, vxc_tau_pw
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      CALL timeset(routineN, handle)

      nspins = SIZE(vxc, 3)
      lsd = (nspins == 2)
      IF (ASSOCIATED(rho_set%rhoa)) THEN
         lsd = .TRUE.
      END IF
      my_fac_triplet = 1.0_dp
      IF (PRESENT(do_triplet)) THEN
         IF (do_triplet) my_fac_triplet = -1.0_dp
      END IF

      CALL xc_rho_set_get(rho_set, drho_cutoff=drho_cutoff)
      xc_fun_section => section_vals_get_subs_vals(xc_section, &
                                                   "XC_FUNCTIONAL")

      !  Calculate the derivatives
      CALL xc_functionals_eval(xc_fun_section, &
                               lsd=lsd, &
                               rho_set=rho_set, &
                               deriv_set=deriv_set, &
                               deriv_order=2)

      CALL divide_by_norm_drho(deriv_set, rho_set, lsd)

      ! multiply by w
      pos => deriv_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         deriv_att%deriv_data(:, :, 1) = w(:, :)*deriv_att%deriv_data(:, :, 1)
      END DO

      NULLIFY (pw_pool)
      ALLOCATE (vxc_pw(nspins))
      DO ispin = 1, nspins
         vxc_pw(ispin)%array => vxc(:, :, ispin:ispin)
      END DO

      NULLIFY (vxc_tau_pw)

      CALL xc_calc_2nd_deriv_analytical(vxc_pw, vxc_tau_pw, deriv_set, rho_set, rho1_set, pw_pool, &
                                        xc_section, gapw=.TRUE., vxg=vxg, tddfpt_fac=my_fac_triplet)

      DEALLOCATE (vxc_pw)

      ! zero the derivative data for the next call
      pos => deriv_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         deriv_att%deriv_data = 0.0_dp
      END DO

      CALL timestop(handle)

   END SUBROUTINE xc_2nd_deriv_of_r

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param needs ...
!> \param nspins ...
!> \param bo ...
! **************************************************************************************************
   SUBROUTINE xc_rho_set_atom_update(rho_set, needs, nspins, bo)

!   This routine allocates the storage arrays for rho and drho
!   In calculate_vxc_atom this is called once for each atomic_kind,
!   After the loop over all the atoms of the kind and over all the points
!   of the radial grid for each atom, rho_set is deallocated.
!   Within the same kind, at each new point on the radial grid, the rho_set
!   arrays rho and drho are overwritten.

      TYPE(xc_rho_set_type), INTENT(INOUT)               :: rho_set
      TYPE(xc_rho_cflags_type), INTENT(IN)               :: needs
      INTEGER, INTENT(IN)                                :: nspins
      INTEGER, DIMENSION(2, 3), INTENT(IN)               :: bo

      INTEGER                                            :: idir

      SELECT CASE (nspins)
      CASE (1)
!     What is this for?
         IF (needs%rho_1_3) THEN
            NULLIFY (rho_set%rho_1_3)
            ALLOCATE (rho_set%rho_1_3(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%rho_1_3 = .TRUE.
            rho_set%has%rho_1_3 = .FALSE.
         END IF
!     Allocate the storage space for the density
         IF (needs%rho) THEN
            NULLIFY (rho_set%rho)
            ALLOCATE (rho_set%rho(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%rho = .TRUE.
            rho_set%has%rho = .FALSE.
         END IF
!     Allocate the storage space for  the norm of the gradient of the density
         IF (needs%norm_drho) THEN
            NULLIFY (rho_set%norm_drho)
            ALLOCATE (rho_set%norm_drho(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%norm_drho = .TRUE.
            rho_set%has%norm_drho = .FALSE.
         END IF
!     Allocate the storage space for the three components of the gradient of the density
         IF (needs%drho) THEN
            DO idir = 1, 3
               NULLIFY (rho_set%drho(idir)%array)
               ALLOCATE (rho_set%drho(idir)%array(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            END DO
            rho_set%owns%drho = .TRUE.
            rho_set%has%drho = .FALSE.
         END IF
      CASE (2)
!     Allocate the storage space for the total density
         IF (needs%rho) THEN
            ! this should never be the case unless you use LDA functionals with LSD
            NULLIFY (rho_set%rho)
            ALLOCATE (rho_set%rho(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%rho = .TRUE.
            rho_set%has%rho = .FALSE.
         END IF
!     What is this for?
         IF (needs%rho_1_3) THEN
            NULLIFY (rho_set%rho_1_3)
            ALLOCATE (rho_set%rho_1_3(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%rho_1_3 = .TRUE.
            rho_set%has%rho_1_3 = .FALSE.
         END IF
!     What is this for?
         IF (needs%rho_spin_1_3) THEN
            NULLIFY (rho_set%rhoa_1_3, rho_set%rhob_1_3)
            ALLOCATE (rho_set%rhoa_1_3(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            ALLOCATE (rho_set%rhob_1_3(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%rho_spin_1_3 = .TRUE.
            rho_set%has%rho_spin_1_3 = .FALSE.
         END IF
!     Allocate the storage space for the spin densities rhoa and rhob
         IF (needs%rho_spin) THEN
            NULLIFY (rho_set%rhoa, rho_set%rhob)
            ALLOCATE (rho_set%rhoa(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            ALLOCATE (rho_set%rhob(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%rho_spin = .TRUE.
            rho_set%has%rho_spin = .FALSE.
         END IF
!     Allocate the storage space for the norm of the gradient of the total density
         IF (needs%norm_drho) THEN
            NULLIFY (rho_set%norm_drho)
            ALLOCATE (rho_set%norm_drho(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%norm_drho = .TRUE.
            rho_set%has%norm_drho = .FALSE.
         END IF
!     Allocate the storage space for the norm of the gradient of rhoa and of rhob separatedly
         IF (needs%norm_drho_spin) THEN
            NULLIFY (rho_set%norm_drhoa, rho_set%norm_drhob)
            ALLOCATE (rho_set%norm_drhoa(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            ALLOCATE (rho_set%norm_drhob(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            rho_set%owns%norm_drho_spin = .TRUE.
            rho_set%has%norm_drho_spin = .FALSE.
         END IF
!     Allocate the storage space for the components of the gradient for the total rho
         IF (needs%drho) THEN
            DO idir = 1, 3
               NULLIFY (rho_set%drho(idir)%array)
               ALLOCATE (rho_set%drho(idir)%array(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            END DO
            rho_set%owns%drho = .TRUE.
            rho_set%has%drho = .FALSE.
         END IF
!     Allocate the storage space for the components of the gradient for rhoa and rhob
         IF (needs%drho_spin) THEN
            DO idir = 1, 3
               NULLIFY (rho_set%drhoa(idir)%array, rho_set%drhob(idir)%array)
               ALLOCATE (rho_set%drhoa(idir)%array(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
               ALLOCATE (rho_set%drhob(idir)%array(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
            END DO
            rho_set%owns%drho_spin = .TRUE.
            rho_set%has%drho_spin = .FALSE.
         END IF
!
      END SELECT

      ! tau part
      IF (needs%tau) THEN
         NULLIFY (rho_set%tau)
         ALLOCATE (rho_set%tau(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         rho_set%owns%tau = .TRUE.
      END IF
      IF (needs%tau_spin) THEN
         NULLIFY (rho_set%tau_a, rho_set%tau_b)
         ALLOCATE (rho_set%tau_a(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         ALLOCATE (rho_set%tau_b(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         rho_set%owns%tau_spin = .TRUE.
         rho_set%has%tau_spin = .FALSE.
      END IF

      ! Laplace part
      IF (needs%laplace_rho) THEN
         NULLIFY (rho_set%laplace_rho)
         ALLOCATE (rho_set%laplace_rho(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         rho_set%owns%laplace_rho = .TRUE.
      END IF
      IF (needs%laplace_rho_spin) THEN
         NULLIFY (rho_set%laplace_rhoa)
         NULLIFY (rho_set%laplace_rhob)
         ALLOCATE (rho_set%laplace_rhoa(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         ALLOCATE (rho_set%laplace_rhob(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         rho_set%owns%laplace_rho_spin = .TRUE.
         rho_set%has%laplace_rho_spin = .TRUE.
      END IF

   END SUBROUTINE xc_rho_set_atom_update

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param lsd ...
!> \param nspins ...
!> \param needs ...
!> \param rho ...
!> \param drho ...
!> \param tau ...
!> \param na ...
!> \param ir ...
! **************************************************************************************************
   SUBROUTINE fill_rho_set(rho_set, lsd, nspins, needs, rho, drho, tau, na, ir)

      TYPE(xc_rho_set_type), INTENT(INOUT)               :: rho_set
      LOGICAL, INTENT(IN)                                :: lsd
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(xc_rho_cflags_type), INTENT(IN)               :: needs
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: rho
      REAL(dp), DIMENSION(:, :, :, :), INTENT(IN)        :: drho
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: tau
      INTEGER, INTENT(IN)                                :: na, ir

      REAL(KIND=dp), PARAMETER                           :: f13 = (1.0_dp/3.0_dp)

      INTEGER                                            :: ia, idir, my_nspins
      LOGICAL                                            :: gradient_f, tddft_split

      my_nspins = nspins
      tddft_split = .FALSE.
      IF (lsd .AND. nspins == 1) THEN
         my_nspins = 2
         tddft_split = .TRUE.
      END IF

      ! some checks
      IF (lsd) THEN
      ELSE
         CPASSERT(SIZE(rho, 3) == 1)
      END IF
      SELECT CASE (my_nspins)
      CASE (1)
         CPASSERT(.NOT. needs%rho_spin)
         CPASSERT(.NOT. needs%drho_spin)
         CPASSERT(.NOT. needs%norm_drho_spin)
         CPASSERT(.NOT. needs%rho_spin_1_3)
      CASE (2)
      CASE default
         CPABORT("Unsupported number of spins")
      END SELECT

      gradient_f = (needs%drho_spin .OR. needs%norm_drho_spin .OR. &
                    needs%drho .OR. needs%norm_drho)

      SELECT CASE (my_nspins)
      CASE (1)
         ! Give rho to 1/3
         IF (needs%rho_1_3) THEN
            DO ia = 1, na
               rho_set%rho_1_3(ia, ir, 1) = MAX(rho(ia, ir, 1), 0.0_dp)**f13
            END DO
            rho_set%owns%rho_1_3 = .TRUE.
            rho_set%has%rho_1_3 = .TRUE.
         END IF
         ! Give the density
         IF (needs%rho) THEN
            DO ia = 1, na
               rho_set%rho(ia, ir, 1) = rho(ia, ir, 1)
            END DO
            rho_set%owns%rho = .TRUE.
            rho_set%has%rho = .TRUE.
         END IF
         ! Give the norm of the gradient of the density
         IF (needs%norm_drho) THEN
            DO ia = 1, na
               rho_set%norm_drho(ia, ir, 1) = drho(4, ia, ir, 1)
            END DO
            rho_set%owns%norm_drho = .TRUE.
            rho_set%has%norm_drho = .TRUE.
         END IF
         ! Give the three components of the gradient of the density
         IF (needs%drho) THEN
            DO idir = 1, 3
               DO ia = 1, na
                  rho_set%drho(idir)%array(ia, ir, 1) = drho(idir, ia, ir, 1)
               END DO
            END DO
            rho_set%owns%drho = .TRUE.
            rho_set%has%drho = .TRUE.
         END IF
      CASE (2)
         ! Give the total density
         IF (needs%rho) THEN
            ! this should never be the case unless you use LDA functionals with LSD
            IF (.NOT. tddft_split) THEN
               DO ia = 1, na
                  rho_set%rho(ia, ir, 1) = rho(ia, ir, 1) + rho(ia, ir, 2)
               END DO
            ELSE
               DO ia = 1, na
                  rho_set%rho(ia, ir, 1) = rho(ia, ir, 1)
               END DO
            END IF
            rho_set%owns%rho = .TRUE.
            rho_set%has%rho = .TRUE.
         END IF
         ! Give the total density to 1/3
         IF (needs%rho_1_3) THEN
            IF (.NOT. tddft_split) THEN
               DO ia = 1, na
                  rho_set%rho_1_3(ia, ir, 1) = MAX(rho(ia, ir, 1) + rho(ia, ir, 2), 0.0_dp)**f13
               END DO
            ELSE
               DO ia = 1, na
                  rho_set%rho_1_3(ia, ir, 1) = MAX(rho(ia, ir, 1), 0.0_dp)**f13
               END DO
            END IF
            rho_set%owns%rho_1_3 = .TRUE.
            rho_set%has%rho_1_3 = .TRUE.
         END IF
         ! Give the spin densities to 1/3
         IF (needs%rho_spin_1_3) THEN
            IF (.NOT. tddft_split) THEN
               DO ia = 1, na
                  rho_set%rhoa_1_3(ia, ir, 1) = MAX(rho(ia, ir, 1), 0.0_dp)**f13
                  rho_set%rhob_1_3(ia, ir, 1) = MAX(rho(ia, ir, 2), 0.0_dp)**f13
               END DO
            ELSE
               DO ia = 1, na
                  rho_set%rhoa_1_3(ia, ir, 1) = MAX(0.5_dp*rho(ia, ir, 1), 0.0_dp)**f13
                  rho_set%rhob_1_3(ia, ir, 1) = rho_set%rhoa_1_3(ia, ir, 1)
               END DO
            END IF
            rho_set%owns%rho_spin_1_3 = .TRUE.
            rho_set%has%rho_spin_1_3 = .TRUE.
         END IF
         ! Give the spin densities rhoa and rhob
         IF (needs%rho_spin) THEN
            IF (.NOT. tddft_split) THEN
               DO ia = 1, na
                  rho_set%rhoa(ia, ir, 1) = rho(ia, ir, 1)
                  rho_set%rhob(ia, ir, 1) = rho(ia, ir, 2)
               END DO
            ELSE
               DO ia = 1, na
                  rho_set%rhoa(ia, ir, 1) = 0.5_dp*rho(ia, ir, 1)
                  rho_set%rhob(ia, ir, 1) = rho_set%rhoa(ia, ir, 1)
               END DO
            END IF
            rho_set%owns%rho_spin = .TRUE.
            rho_set%has%rho_spin = .TRUE.
         END IF
         ! Give the norm of the gradient of the total density
         IF (needs%norm_drho) THEN
            IF (.NOT. tddft_split) THEN
               DO ia = 1, na
                  rho_set%norm_drho(ia, ir, 1) = SQRT( &
                                                 (drho(1, ia, ir, 1) + drho(1, ia, ir, 2))**2 + &
                                                 (drho(2, ia, ir, 1) + drho(2, ia, ir, 2))**2 + &
                                                 (drho(3, ia, ir, 1) + drho(3, ia, ir, 2))**2)
               END DO
            ELSE
               DO ia = 1, na
                  rho_set%norm_drho(ia, ir, 1) = drho(4, ia, ir, 1)
               END DO
            END IF
            rho_set%owns%norm_drho = .TRUE.
            rho_set%has%norm_drho = .TRUE.
         END IF
         ! Give the norm of the gradient of rhoa and of rhob separatedly
         IF (needs%norm_drho_spin) THEN
            IF (.NOT. tddft_split) THEN
               DO ia = 1, na
                  rho_set%norm_drhoa(ia, ir, 1) = drho(4, ia, ir, 1)
                  rho_set%norm_drhob(ia, ir, 1) = drho(4, ia, ir, 2)
               END DO
            ELSE
               DO ia = 1, na
                  rho_set%norm_drhoa(ia, ir, 1) = 0.5_dp*drho(4, ia, ir, 1)
                  rho_set%norm_drhob(ia, ir, 1) = rho_set%norm_drhoa(ia, ir, 1)
               END DO
            END IF
            rho_set%owns%norm_drho_spin = .TRUE.
            rho_set%has%norm_drho_spin = .TRUE.
         END IF
         ! Give the components of the gradient for the total rho
         IF (needs%drho) THEN
            IF (.NOT. tddft_split) THEN
               DO idir = 1, 3
                  DO ia = 1, na
                     rho_set%drho(idir)%array(ia, ir, 1) = drho(idir, ia, ir, 1) + drho(idir, ia, ir, 2)
                  END DO
               END DO
            ELSE
               DO idir = 1, 3
                  DO ia = 1, na
                     rho_set%drho(idir)%array(ia, ir, 1) = drho(idir, ia, ir, 1)
                  END DO
               END DO
            END IF
            rho_set%owns%drho = .TRUE.
            rho_set%has%drho = .TRUE.
         END IF
         ! Give the components of the gradient for rhoa and rhob
         IF (needs%drho_spin) THEN
            IF (.NOT. tddft_split) THEN
               DO idir = 1, 3
                  DO ia = 1, na
                     rho_set%drhoa(idir)%array(ia, ir, 1) = drho(idir, ia, ir, 1)
                     rho_set%drhob(idir)%array(ia, ir, 1) = drho(idir, ia, ir, 2)
                  END DO
               END DO
            ELSE
               DO idir = 1, 3
                  DO ia = 1, na
                     rho_set%drhoa(idir)%array(ia, ir, 1) = 0.5_dp*drho(idir, ia, ir, 1)
                     rho_set%drhob(idir)%array(ia, ir, 1) = rho_set%drhoa(idir)%array(ia, ir, 1)
                  END DO
               END DO
            END IF
            rho_set%owns%drho_spin = .TRUE.
            rho_set%has%drho_spin = .TRUE.
         END IF
         !
      END SELECT

      ! tau part
      IF (needs%tau .OR. needs%tau_spin) THEN
         CPASSERT(SIZE(tau, 3) == my_nspins)
      END IF
      IF (needs%tau) THEN
         IF (my_nspins == 2) THEN
            DO ia = 1, na
               rho_set%tau(ia, ir, 1) = tau(ia, ir, 1) + tau(ia, ir, 2)
            END DO
            rho_set%owns%tau = .TRUE.
            rho_set%has%tau = .TRUE.
         ELSE
            DO ia = 1, na
               rho_set%tau(ia, ir, 1) = tau(ia, ir, 1)
            END DO
            rho_set%owns%tau = .TRUE.
            rho_set%has%tau = .TRUE.
         END IF
      END IF
      IF (needs%tau_spin) THEN
         DO ia = 1, na
            rho_set%tau_a(ia, ir, 1) = tau(ia, ir, 1)
            rho_set%tau_b(ia, ir, 1) = tau(ia, ir, 2)
         END DO
         rho_set%owns%tau_spin = .TRUE.
         rho_set%has%tau_spin = .TRUE.
      END IF

   END SUBROUTINE fill_rho_set

END MODULE xc_atom
