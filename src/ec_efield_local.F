!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates the energy contribution and the mo_derivative of
!>        a static electric field (nonperiodic)
!> \par History
!>      Adjusted from qs_efield_local
!> \author JGH (10.2019)
! **************************************************************************************************
MODULE ec_efield_local
   USE ai_moments,                      ONLY: dipole_force
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE ec_env_types,                    ONLY: energy_correction_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_period_efield_types,          ONLY: efield_berry_type,&
                                              init_efield_matrices,&
                                              set_efield_matrices
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_efield_local'

   ! *** Public subroutines ***

   PUBLIC :: ec_efield_local_operator, ec_efield_integrals

! **************************************************************************************************

CONTAINS

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE ec_efield_local_operator(qs_env, ec_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_efield_local_operator'

      INTEGER                                            :: handle
      REAL(dp), DIMENSION(3)                             :: rpoint
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (dft_control%apply_efield) THEN
         rpoint = 0.0_dp
         CALL ec_efield_integrals(qs_env, ec_env, rpoint)
         CALL ec_efield_mo_derivatives(qs_env, ec_env, rpoint, calculate_forces)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_efield_local_operator

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param rpoint ...
! **************************************************************************************************
   SUBROUTINE ec_efield_integrals(qs_env, ec_env, rpoint)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rpoint

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_efield_integrals'

      INTEGER                                            :: handle, i
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_s
      TYPE(efield_berry_type), POINTER                   :: efield, efieldref

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, efield=efieldref)
      efield => ec_env%efield
      CALL init_efield_matrices(efield)
      matrix_s => ec_env%matrix_s(:, 1)
      ALLOCATE (dipmat(3))
      DO i = 1, 3
         ALLOCATE (dipmat(i)%matrix)
         CALL dbcsr_copy(dipmat(i)%matrix, matrix_s(1)%matrix, 'DIP MAT')
         CALL dbcsr_set(dipmat(i)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, dipmat, 1, rpoint, basis_type="HARRIS")
      CALL set_efield_matrices(efield=efield, dipmat=dipmat)
      ec_env%efield => efield

      CALL timestop(handle)

   END SUBROUTINE ec_efield_integrals

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param rpoint ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE ec_efield_mo_derivatives(qs_env, ec_env, rpoint, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rpoint
      LOGICAL                                            :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_efield_mo_derivatives'

      INTEGER :: atom_a, atom_b, handle, i, ia, iatom, icol, idir, ikind, irow, iset, ispin, &
         jatom, jkind, jset, ldab, natom, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found, trans
      REAL(dp)                                           :: charge, dab, fdir
      REAL(dp), DIMENSION(3)                             :: ci, fieldpol, ra, rab, rac, rbc, ria
      REAL(dp), DIMENSION(3, 3)                          :: forcea, forceb
      REAL(dp), DIMENSION(:, :), POINTER                 :: p_block_a, p_block_b, pblock, pmat, work
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(efield_berry_type), POINTER                   :: efield
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, cell=cell, particle_set=particle_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      energy=energy, para_env=para_env, sab_orb=sab_orb)

      efield => ec_env%efield

      fieldpol = dft_control%efield_fields(1)%efield%polarisation* &
                 dft_control%efield_fields(1)%efield%strength

      ! nuclear contribution
      natom = SIZE(particle_set)
      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force)
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
      END IF
      ci = 0.0_dp
      DO ia = 1, natom
         CALL get_atomic_kind(particle_set(ia)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge)
         ria = particle_set(ia)%r - rpoint
         ria = pbc(ria, cell)
         ci(:) = ci(:) + charge*ria(:)
         IF (calculate_forces) THEN
            IF (para_env%mepos == 0) THEN
               iatom = atom_of_kind(ia)
               DO idir = 1, 3
                  force(ikind)%efield(idir, iatom) = force(ikind)%efield(idir, iatom) - fieldpol(idir)*charge
               END DO
            END IF
         END IF
      END DO

      IF (ec_env%should_update) THEN
         ec_env%efield_nuclear = -SUM(ci(:)*fieldpol(:))
         ! Update KS matrix
         matrix_ks => ec_env%matrix_h(:, 1)
         dipmat => efield%dipmat
         DO ispin = 1, SIZE(matrix_ks)
            DO idir = 1, 3
               CALL dbcsr_add(matrix_ks(ispin)%matrix, dipmat(idir)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=fieldpol(idir))
            END DO
         END DO
      END IF

      ! forces from the efield contribution
      IF (calculate_forces) THEN
         nkind = SIZE(qs_kind_set)
         natom = SIZE(particle_set)

         ALLOCATE (basis_set_list(nkind))
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type="HARRIS")
            IF (ASSOCIATED(basis_set_a)) THEN
               basis_set_list(ikind)%gto_basis_set => basis_set_a
            ELSE
               NULLIFY (basis_set_list(ikind)%gto_basis_set)
            END IF
         END DO
         !
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rab)
            basis_set_a => basis_set_list(ikind)%gto_basis_set
            IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
            basis_set_b => basis_set_list(jkind)%gto_basis_set
            IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
            ! basis ikind
            first_sgfa => basis_set_a%first_sgf
            la_max => basis_set_a%lmax
            la_min => basis_set_a%lmin
            npgfa => basis_set_a%npgf
            nseta = basis_set_a%nset
            nsgfa => basis_set_a%nsgf_set
            rpgfa => basis_set_a%pgf_radius
            set_radius_a => basis_set_a%set_radius
            sphi_a => basis_set_a%sphi
            zeta => basis_set_a%zet
            ! basis jkind
            first_sgfb => basis_set_b%first_sgf
            lb_max => basis_set_b%lmax
            lb_min => basis_set_b%lmin
            npgfb => basis_set_b%npgf
            nsetb = basis_set_b%nset
            nsgfb => basis_set_b%nsgf_set
            rpgfb => basis_set_b%pgf_radius
            set_radius_b => basis_set_b%set_radius
            sphi_b => basis_set_b%sphi
            zetb => basis_set_b%zet

            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)

            ra(:) = particle_set(iatom)%r(:) - rpoint(:)
            rac(:) = pbc(ra(:), cell)
            rbc(:) = rac(:) + rab(:)
            dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
               trans = .FALSE.
            ELSE
               irow = jatom
               icol = iatom
               trans = .TRUE.
            END IF

            fdir = 2.0_dp
            IF (iatom == jatom .AND. dab < 1.e-10_dp) fdir = 1.0_dp

            ! density matrix
            NULLIFY (p_block_a)
            CALL dbcsr_get_block_p(ec_env%matrix_p(1, 1)%matrix, irow, icol, p_block_a, found)
            IF (.NOT. found) CYCLE
            IF (SIZE(ec_env%matrix_p, 1) > 1) THEN
               NULLIFY (p_block_b)
               CALL dbcsr_get_block_p(ec_env%matrix_p(2, 1)%matrix, irow, icol, p_block_b, found)
               CPASSERT(found)
            END IF
            forcea = 0.0_dp
            forceb = 0.0_dp

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               DO jset = 1, nsetb
                  IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)
                  ! Calculate the primitive integrals (da|O|b) and (a|O|db)
                  ldab = MAX(ncoa, ncob)
                  ALLOCATE (work(ldab, ldab), pmat(ncoa, ncob))
                  ! Decontract P matrix block
                  pmat = 0.0_dp
                  DO i = 1, SIZE(ec_env%matrix_p, 1)
                     IF (i == 1) THEN
                        pblock => p_block_a
                     ELSE
                        pblock => p_block_b
                     END IF
                     IF (.NOT. ASSOCIATED(pblock)) CYCLE
                     IF (trans) THEN
                        CALL dgemm("N", "T", ncoa, nsgfb(jset), nsgfa(iset), &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   pblock(sgfb, sgfa), SIZE(pblock, 1), &
                                   0.0_dp, work(1, 1), ldab)
                     ELSE
                        CALL dgemm("N", "N", ncoa, nsgfb(jset), nsgfa(iset), &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   pblock(sgfa, sgfb), SIZE(pblock, 1), &
                                   0.0_dp, work(1, 1), ldab)
                     END IF
                     CALL dgemm("N", "T", ncoa, ncob, nsgfb(jset), &
                                1.0_dp, work(1, 1), ldab, &
                                sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                                1.0_dp, pmat(1, 1), ncoa)
                  END DO

                  CALL dipole_force(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                    lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                    1, rac, rbc, pmat, forcea, forceb)

                  DEALLOCATE (work, pmat)
               END DO
            END DO

            DO idir = 1, 3
               force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) &
                                                  + fdir*fieldpol(idir)*forcea(idir, 1:3)
               force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) &
                                                  + fdir*fieldpol(idir)*forceb(idir, 1:3)
            END DO

         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
         DEALLOCATE (basis_set_list)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_efield_mo_derivatives

END MODULE ec_efield_local
