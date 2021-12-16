!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Simplified Tamm Dancoff approach (sTDA).
! **************************************************************************************************
MODULE qs_tddfpt2_stda_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE basis_set_container_types,       ONLY: get_basis_from_container
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              stda_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_row_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_schur_product
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_element, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_to_fm_submat, &
        cp_fm_type, cp_fm_vectorssum
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iterate
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_create, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE ewald_environment_types,         ONLY: ewald_env_create,&
                                              ewald_env_get,&
                                              ewald_env_set,&
                                              ewald_environment_type,&
                                              read_ewald_section_tb
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_evaluate
   USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                              ewald_pw_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos,&
                                              tddfpt_work_matrices
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_stda_utils'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC:: stda_init_matrices, stda_calculate_kernel, get_lowdin_x, setup_gamma

CONTAINS

! **************************************************************************************************
!> \brief Calculate sTDA matrices
!> \param qs_env ...
!> \param stda_kernel ...
!> \param sub_env ...
!> \param work ...
!> \param tddfpt_control ...
! **************************************************************************************************
   SUBROUTINE stda_init_matrices(qs_env, stda_kernel, sub_env, work, tddfpt_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type)                                :: stda_kernel
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_init_matrices'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_coulomb
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section

      CALL timeset(routineN, handle)

      do_coulomb = .NOT. tddfpt_control%rks_triplets
      IF (do_coulomb) THEN
         ! calculate exchange gamma matrix
         CALL setup_gamma(qs_env, stda_kernel, sub_env, work%gamma_exchange)
      END IF

      ! calculate S_half and Lowdin MO coefficients
      CALL get_lowdin_mo_coefficients(qs_env, sub_env, work)

      ! initialize Ewald for sTDA
      IF (tddfpt_control%stda_control%do_ewald) THEN
         NULLIFY (ewald_env, ewald_pw)
         ALLOCATE (ewald_env)
         CALL ewald_env_create(ewald_env, sub_env%para_env)
         poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
         CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
         ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
         print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
         CALL get_qs_env(qs_env, cell=cell, cell_ref=cell_ref)
         CALL read_ewald_section_tb(ewald_env, ewald_section, cell_ref%hmat)
         ALLOCATE (ewald_pw)
         CALL ewald_pw_create(ewald_pw, ewald_env, cell, cell_ref, print_section=print_section)
         work%ewald_env => ewald_env
         work%ewald_pw => ewald_pw
      END IF

      CALL timestop(handle)

   END SUBROUTINE stda_init_matrices
! **************************************************************************************************
!> \brief Calculate sTDA exchange-type contributions
!> \param qs_env ...
!> \param stda_env ...
!> \param sub_env ...
!> \param gamma_matrix sTDA exchange-type contributions
!> \param ndim ...
!> \note  Note the specific sTDA notation exchange-type integrals (ia|jb) refer to Coulomb interaction
! **************************************************************************************************
   SUBROUTINE setup_gamma(qs_env, stda_env, sub_env, gamma_matrix, ndim)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type)                                :: stda_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: ndim

      CHARACTER(len=*), PARAMETER                        :: routineN = 'setup_gamma'
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      INTEGER                                            :: handle, i, iatom, icol, ikind, imat, &
                                                            irow, jatom, jkind, natom, nmat
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dfcut, dgb, dr, eta, fcut, r, rcut, &
                                                            rcuta, rcutb, x
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dgblock, gblock
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, natom=natom)
      dbcsr_dist => sub_env%dbcsr_dist
      ! Using the overlap list here can have a considerable effect on the number of
      ! terms calculated. This makes gamma also dependent on EPS_DEFAULT -> Overlap
      n_list => sub_env%sab_orb

      IF (PRESENT(ndim)) THEN
         nmat = ndim
      ELSE
         nmat = 1
      END IF
      CPASSERT(nmat == 1 .OR. nmat == 4)
      CPASSERT(.NOT. ASSOCIATED(gamma_matrix))
      CALL dbcsr_allocate_matrix_set(gamma_matrix, nmat)

      ALLOCATE (row_blk_sizes(natom))
      row_blk_sizes(1:natom) = 1
      DO imat = 1, nmat
         ALLOCATE (gamma_matrix(imat)%matrix)
      END DO

      CALL dbcsr_create(gamma_matrix(1)%matrix, name="gamma", dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_symmetric, row_blk_size=row_blk_sizes, &
                        col_blk_size=row_blk_sizes, nze=0)
      DO imat = 2, nmat
         CALL dbcsr_create(gamma_matrix(imat)%matrix, name="dgamma", dist=dbcsr_dist, &
                           matrix_type=dbcsr_type_antisymmetric, row_blk_size=row_blk_sizes, &
                           col_blk_size=row_blk_sizes, nze=0)
      END DO

      DEALLOCATE (row_blk_sizes)

      ! setup the matrices using the neighbor list
      DO imat = 1, nmat
         CALL cp_dbcsr_alloc_block_from_nbl(gamma_matrix(imat)%matrix, n_list)
         CALL dbcsr_set(gamma_matrix(imat)%matrix, 0.0_dp)
      END DO

      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         dr = SQRT(SUM(rij(:)**2)) ! interatomic distance

         eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param + &
                stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         NULLIFY (gblock)
         CALL dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix, &
                                row=irow, col=icol, BLOCK=gblock, found=found)
         CPASSERT(found)

         ! get rcuta and rcutb
         rcuta = stda_env%kind_param_set(ikind)%kind_param%rcut
         rcutb = stda_env%kind_param_set(jkind)%kind_param%rcut
         rcut = rcuta + rcutb

         !>   Computes the short-range gamma parameter from
         !>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
         IF (dr < 1.e-6) THEN
            ! on site terms
            gblock(:, :) = gblock(:, :) + eta
         ELSEIF (dr > rcut) THEN
            ! do nothing
         ELSE
            IF (dr < rcut - rsmooth) THEN
               fcut = 1.0_dp
            ELSE
               r = dr - (rcut - rsmooth)
               x = r/rsmooth
               fcut = -6._dp*x**5 + 15._dp*x**4 - 10._dp*x**3 + 1._dp
            END IF
            gblock(:, :) = gblock(:, :) + &
                           fcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                           **(1._dp/stda_env%alpha_param) - fcut/dr
         END IF

         IF (nmat > 1) THEN
            !>   Computes the short-range gamma parameter from
            !>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
            !>   Derivatives
            IF (dr < 1.e-6 .OR. dr > rcut) THEN
               ! on site terms or beyond cutoff
               dgb = 0.0_dp
            ELSE
               IF (dr < rcut - rsmooth) THEN
                  fcut = 1.0_dp
                  dfcut = 0.0_dp
               ELSE
                  r = dr - (rcut - rsmooth)
                  x = r/rsmooth
                  fcut = -6._dp*x**5 + 15._dp*x**4 - 10._dp*x**3 + 1._dp
                  dfcut = -30._dp*x**4 + 60._dp*x**3 - 30._dp*x**2
                  dfcut = dfcut/rsmooth
               END IF
               dgb = dfcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                     **(1._dp/stda_env%alpha_param)
               dgb = dgb - dfcut/dr + fcut/dr**2
               dgb = dgb - fcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                     **(1._dp/stda_env%alpha_param + 1._dp)*dr**(stda_env%alpha_param - 1._dp)
            END IF
            DO imat = 2, nmat
               NULLIFY (dgblock)
               CALL dbcsr_get_block_p(matrix=gamma_matrix(imat)%matrix, &
                                      row=irow, col=icol, BLOCK=dgblock, found=found)
               IF (found) THEN
                  IF (dr > 1.e-6) THEN
                     i = imat - 1
                     IF (irow == iatom) THEN
                        dgblock(:, :) = dgblock(:, :) + dgb*rij(i)/dr
                     ELSE
                        dgblock(:, :) = dgblock(:, :) - dgb*rij(i)/dr
                     END IF
                  END IF
               END IF
            END DO
         END IF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      DO imat = 1, nmat
         CALL dbcsr_finalize(gamma_matrix(imat)%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE setup_gamma

! **************************************************************************************************
!> \brief Calculate Lowdin MO coefficients
!> \param qs_env ...
!> \param sub_env ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_mo_coefficients(qs_env, sub_env, work)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_mo_coefficients'

      INTEGER                                            :: handle, i, iounit, ispin, j, &
                                                            max_iter_lanczos, nactive, ndep, nsgf, &
                                                            nspins, order_lanczos
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: eps_lanczos, sij, threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: slam
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type)                                   :: fm_s_half, fm_work1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type)                                   :: sm_hinv
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_s
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      ! Calculate S^1/2 matrix
      IF (iounit > 0) THEN
         WRITE (iounit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                      sTDA Kernel: Create Matrix SQRT(S)                     -", &
            "-------------------------------------------------------------------------------"
      END IF

      IF (sub_env%is_split) THEN
         CPABORT('SPLIT')
      ELSE
         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrixkp_s)
         CPASSERT(ASSOCIATED(matrixkp_s))
         IF (SIZE(matrixkp_s, 2) > 1) CPWARN("not implemented for k-points.")
         sm_s => matrixkp_s(1, 1)%matrix
      END IF
      sm_h => work%shalf

      CALL dbcsr_create(sm_hinv, template=sm_s)
      CALL dbcsr_add_on_diag(sm_h, 1.0_dp)
      threshold = 1.0e-8_dp
      order_lanczos = 3
      eps_lanczos = 1.0e-4_dp
      max_iter_lanczos = 40
      CALL matrix_sqrt_Newton_Schulz(sm_h, sm_hinv, sm_s, &
                                     threshold, order_lanczos, eps_lanczos, max_iter_lanczos, &
                                     converged=converged)
      CALL dbcsr_release(sm_hinv)
      !
      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      ! Get the total number of contracted spherical Gaussian basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
      !
      IF (.NOT. converged) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, "(T3,A)") "STDA| Newton-Schulz iteration did not converge"
            WRITE (iounit, "(T3,A)") "STDA| Calculate SQRT(S) from diagonalization"
         END IF
         CALL get_qs_env(qs_env=qs_env, scf_control=scf_control)
         ! Provide full size work matrices
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=sub_env%para_env, &
                                  context=sub_env%blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=nsgf)
         CALL cp_fm_create(matrix=fm_s_half, matrix_struct=fmstruct, name="S^(1/2) MATRIX")
         CALL cp_fm_create(matrix=fm_work1, matrix_struct=fmstruct, name="TMP MATRIX")
         CALL cp_fm_struct_release(fmstruct=fmstruct)
         CALL copy_dbcsr_to_fm(sm_s, fm_s_half)
         CALL cp_fm_power(fm_s_half, fm_work1, 0.5_dp, scf_control%eps_eigval, ndep)
         IF (ndep /= 0) &
            CALL cp_warn(__LOCATION__, &
                         "Overlap matrix exhibits linear dependencies. At least some "// &
                         "eigenvalues have been quenched.")
         CALL copy_fm_to_dbcsr(fm_s_half, sm_h)
         CALL cp_fm_release(matrix=fm_s_half)
         CALL cp_fm_release(matrix=fm_work1)
         IF (iounit > 0) WRITE (iounit, *)
      END IF

      nspins = SIZE(sub_env%mos_occ)

      DO ispin = 1, nspins
         CALL cp_fm_get_info(work%ctransformed(ispin)%matrix, ncol_global=nactive)
         CALL cp_dbcsr_sm_fm_multiply(work%shalf, sub_env%mos_occ(ispin)%matrix, &
                                      work%ctransformed(ispin)%matrix, nactive, alpha=1.0_dp, beta=0.0_dp)
      END DO

      ! for Lowdin forces
      CALL cp_fm_create(matrix=fm_work1, matrix_struct=work%S_eigenvectors%matrix_struct, name="TMP MATRIX")
      CALL copy_dbcsr_to_fm(sm_s, fm_work1)
      CALL choose_eigv_solver(fm_work1, work%S_eigenvectors, work%S_eigenvalues)
      CALL cp_fm_release(matrix=fm_work1)
      !
      ALLOCATE (slam(nsgf, 1))
      DO i = 1, nsgf
         IF (work%S_eigenvalues(i) > 0._dp) THEN
            slam(i, 1) = SQRT(work%S_eigenvalues(i))
         ELSE
            CPABORT("S matrix not positive definit")
         END IF
      END DO
      DO i = 1, nsgf
         CALL cp_fm_set_submatrix(work%slambda, slam, 1, i, nsgf, 1, 1.0_dp, 0.0_dp)
      END DO
      DO i = 1, nsgf
         CALL cp_fm_set_submatrix(work%slambda, slam, i, 1, 1, nsgf, 1.0_dp, 1.0_dp, .TRUE.)
      END DO
      CALL cp_fm_get_info(work%slambda, local_data=local_data)
      DO i = 1, SIZE(local_data, 2)
         DO j = 1, SIZE(local_data, 1)
            sij = local_data(j, i)
            IF (sij > 0.0_dp) sij = 1.0_dp/sij
            local_data(j, i) = sij
         END DO
      END DO
      DEALLOCATE (slam)

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_mo_coefficients

! **************************************************************************************************
!> \brief Calculate Lowdin transformed Davidson trial vector X
!>        shalf (dbcsr), xvec, xt (fm) are defined in the same sub_env
!> \param shalf ...
!> \param xvec ...
!> \param xt ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_x(shalf, xvec, xt)

      TYPE(dbcsr_type)                                   :: shalf
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: xvec
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: xt

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_lowdin_x'

      INTEGER                                            :: handle, ispin, nactive, nspins

      CALL timeset(routineN, handle)

      nspins = SIZE(xvec)

      ! Build Lowdin transformed tilde(X)= S^1/2 X for each spin
      DO ispin = 1, nspins
         CALL cp_fm_get_info(xt(ispin)%matrix, ncol_global=nactive)
         CALL cp_dbcsr_sm_fm_multiply(shalf, xvec(ispin)%matrix, &
                                      xt(ispin)%matrix, nactive, alpha=1.0_dp, beta=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_x
! **************************************************************************************************
!> \brief More output for sTDA when switching on keyword SHIFT_OUTPUT
!> \brief Prints gap
!> \param qs_env ...
!> \param sub_env ...
!> \param stda_control ...
!> \param work ...
!> \param gs_mos ...
!> \param ivect ...
! **************************************************************************************************
   SUBROUTINE more_output(qs_env, sub_env, stda_control, work, gs_mos, ivect)
!
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(stda_control_type)                            :: stda_control
      TYPE(tddfpt_work_matrices), INTENT(IN)             :: work
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: ivect

      CHARACTER(len=*), PARAMETER                        :: routineN = 'more_output'

      LOGICAL                                            :: lshift_output
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

!
!!

      INTEGER:: ispin, start, handle, ii, iounit, n_atoms, nspins, limit, nao, itertddft
      INTEGER, DIMENSION(maxspins)                       :: nactive, nvir, nmo

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=itertddft)

      lshift_output = stda_control%shift_output

      IF ((itertddft .EQ. 1) .AND. (ivect .EQ. 1)) THEN

         nspins = SIZE(gs_mos)
         DO ispin = 1, nspins
            CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, ncol_global=nvir(ispin), nrow_global=nao)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nactive(ispin))
            nmo(ispin) = nactive(ispin) + nvir(ispin)
         END DO

         NULLIFY (particle_set)
         CALL get_qs_env(qs_env=qs_env, &
                         particle_set=particle_set)
         para_env => sub_env%para_env
         blacs_env => sub_env%blacs_env
         n_atoms = SIZE(particle_set)

         IF (iounit > 0) WRITE (iounit, '(/,1X,A,/)') "Frontier orbitals in a.u. / eV"
         DO ispin = 1, nspins
            IF ((iounit > 0) .AND. (nspins .EQ. 2)) WRITE (iounit, '(A,I4)') "Spin", ispin
            limit = MIN(nactive(ispin), 10)
            start = nactive(ispin) - limit + 1
            DO ii = start, nactive(ispin)
               IF (ii .EQ. nactive(ispin)) THEN
                  IF (iounit > 0) WRITE (iounit, "(1X,A,3X,I4,3X,F10.8,3X,F15.6)") "HOMO:", ii, gs_mos(ispin)%evals_occ(ii), &
                     gs_mos(ispin)%evals_occ(ii)*evolt
               ELSE
                  IF (iounit > 0) WRITE (iounit, "(9X,I4,3X,F10.8,3X,F15.6)") ii, gs_mos(ispin)%evals_occ(ii), &
                     gs_mos(ispin)%evals_occ(ii)*evolt
               END IF
            END DO
            limit = MIN(nvir(ispin), 10)
            DO ii = 1, limit
               IF (ii .EQ. 1) THEN
      IF (iounit > 0) WRITE (iounit, "(1X,A,3X,I4,3X,F10.8,3X,F15.6)") "LUMO:", nactive(ispin) + ii, gs_mos(ispin)%evals_virt(ii), &
                     gs_mos(ispin)%evals_virt(ii)*evolt
               ELSE
                 IF (iounit > 0) WRITE (iounit, "(1X,8X,I4,3X,F10.8,3X,F15.6)") nactive(ispin) + ii, gs_mos(ispin)%evals_virt(ii), &
                     gs_mos(ispin)%evals_virt(ii)*evolt
               END IF
            END DO

            IF (iounit > 0) WRITE (iounit, '(/,1X,A,3X,F15.6,3X,F15.6)') "HOMO-LUMO gap [ a.u. / eV ]", &
               gs_mos(ispin)%evals_virt(1) - gs_mos(ispin)%evals_occ(nactive(ispin)), &
               (gs_mos(ispin)%evals_virt(1) - gs_mos(ispin)%evals_occ(nactive(ispin)))*evolt
            IF (iounit > 0) WRITE (iounit, '(/,1X,A,3X,F15.6,3X,F15.6,/,/)') "Fermi level [ a.u. / eV ]", &
               gs_mos(ispin)%evals_occ(nactive(ispin)), gs_mos(ispin)%evals_occ(nactive(ispin))*evolt

         END DO

         DO ispin = 1, nspins
            NULLIFY (fmstruct)
            CALL cp_fm_get_info(work%ctransformed(ispin)%matrix, matrix_struct=fmstruct)
         END DO

      END IF !itertddft

      CALL timestop(handle)
!
   END SUBROUTINE
! **************************************************************************************************
!> \brief Calculate shift for diagonal kernel contributions for xTB-sTDA as suggested originally
!>        by Grimme
!> \param qs_env ...
!> \param sub_env ...
!> \param stda_control ...
!> \param work ...
!> \param gs_mos ...
!> \param X ...
!> \param ivect ...
!> \param res ...
! **************************************************************************************************
   SUBROUTINE k_shift_orig(qs_env, sub_env, stda_control, work, gs_mos, X, ivect, res)
!
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(stda_control_type)                            :: stda_control
      TYPE(tddfpt_work_matrices), INTENT(IN)             :: work
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: X
      INTEGER, INTENT(in)                                :: ivect
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: res

      CHARACTER(len=*), PARAMETER                        :: routineN = 'k_shift_orig'

      INTEGER                                            :: aa, blk, dimen, ewald_type, handle, &
                                                            iatom, ii, ikind, iocc, iounit, ispin, &
                                                            itertddft, ivir, jatom, jj, kk, &
                                                            n_atoms, nao, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(2)                              :: nlim
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      INTEGER, DIMENSION(maxspins)                       :: nactive, nmo, nvir
      LOGICAL                                            :: calculate_forces, ldo_ewald, &
                                                            lshift_output, use_virial
      REAL(KIND=dp) :: alpha, dr, factor, factortwo, prefactor, shift_delta, shift_os_delta, &
         shift_os_sigma, shift_sigma, spinfac, steepness, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: charge_per_atom, csquared
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: charge_per_atom_2index, &
                                                            gcharge_per_atom, &
                                                            gcharge_per_atom_2index
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: ctransformed_virtual, DeltaK, diagshift, &
                                                            diagshiftmo, SC, Xmo
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: charges, gcharges
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstructia
      TYPE(cp_fm_type), POINTER                          :: cmatocc, cmatvir, tmpmat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type), POINTER                          :: sm_s, tempmat
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_iatom
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(virial_type), POINTER                         :: virial

!
!

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=itertddft)

      shift_delta = stda_control%diag_shift_delta/evolt
      shift_sigma = stda_control%diag_shift_sigma/evolt
      shift_os_delta = stda_control%diag_os_shift_delta/evolt
      shift_os_sigma = stda_control%diag_os_shift_sigma/evolt

      ldo_ewald = stda_control%do_ewald
      lshift_output = stda_control%shift_output

      NULLIFY (particle_set, atomic_kind_set, kind_of, qs_kind_set, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell)

      n_atoms = SIZE(particle_set)

      ALLOCATE (kind_of(n_atoms))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      ALLOCATE (first_sgf(n_atoms))
      ALLOCATE (last_sgf(n_atoms))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      para_env => sub_env%para_env
      blacs_env => sub_env%blacs_env

      nspins = SIZE(gs_mos)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, ncol_global=nvir(ispin), nrow_global=nao)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nactive(ispin))
         nmo(ispin) = nactive(ispin) + nvir(ispin)
      END DO

      spinfac = 2.0_dp
      IF (nspins == 2) spinfac = 1.0_dp
!
! calculate Ctransformed_virtual
!
      NULLIFY (ctransformed_virtual)
      ALLOCATE (ctransformed_virtual(nspins))
      DO ispin = 1, nspins
         NULLIFY (fmstruct)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=ctransformed_virtual(ispin)%matrix, matrix_struct=fmstruct)
         CALL cp_dbcsr_sm_fm_multiply(work%shalf, gs_mos(ispin)%mos_virt, &
                                      ctransformed_virtual(ispin)%matrix, nvir(ispin), alpha=1.0_dp, beta=0.0_dp)
         IF (debug_this_module) THEN
!    debugging test: Calculating Loewdin transition density charges for the complete mo basis
            dimen = nmo(ispin)
            ALLOCATE (csquared(dimen))
            csquared = 0.0_dp
            DO ii = 1, nactive(ispin)
               DO jj = 1, nao
                  CALL cp_fm_get_element(ctransformed_virtual(ispin)%matrix, jj, ii, tmp)
                  csquared(ii) = csquared(ii) + (ABS(tmp))**2
               END DO
            END DO
            DO ii = 1, nvir(ispin)
               DO jj = 1, nao
                  CALL cp_fm_get_element(ctransformed_virtual(ispin)%matrix, jj, ii, tmp)
                  csquared(nactive(ispin) + ii) = csquared(nactive(ispin) + ii) + (ABS(tmp))**2
               END DO
            END DO
            IF (debug_this_module .AND. (iounit > 0)) THEN
               IF ((itertddft .EQ. 1) .AND. (ivect .EQ. 1)) THEN
               DO ii = 1, nmo(ispin)
                  WRITE (UNIT=iounit, FMT="(A,F10.8)") 'check sum loewdin coefficient != 1  :  ', csquared(ii)
               END DO
               END IF
            END IF
            DEALLOCATE (csquared)
         END IF
      END DO !nspins

! get charges  q_ia^A
      NULLIFY (charges)
      ALLOCATE (charges(nspins, n_atoms))
      NULLIFY (gcharges)
      ALLOCATE (gcharges(nspins, n_atoms))
      DO ispin = 1, nspins
         NULLIFY (fmstructia)
         CALL cp_fm_struct_create(fmstruct=fmstructia, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nvir(ispin), &
                                  ncol_global=nactive(ispin))
         DO iatom = 1, n_atoms
            CALL cp_fm_create(matrix=charges(ispin, iatom)%matrix, matrix_struct=fmstructia)
            CALL cp_fm_create(matrix=gcharges(ispin, iatom)%matrix, matrix_struct=fmstructia)
            CALL cp_fm_set_all(gcharges(ispin, iatom)%matrix, 0.0_dp)

            NULLIFY (basis_set_iatom, qs_kind)
            ikind = kind_of(iatom)
            qs_kind => qs_kind_set(ikind)
            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set_iatom, basis_type="ORB")

            NULLIFY (cmatocc, cmatvir, fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set_iatom%nsgf, &
                                     ncol_global=nactive(ispin))
            CALL cp_fm_create(matrix=cmatocc, matrix_struct=fmstruct)
            CALL cp_fm_struct_release(fmstruct=fmstruct)

            NULLIFY (fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set_iatom%nsgf, &
                                     ncol_global=nvir(ispin))
            CALL cp_fm_create(matrix=cmatvir, matrix_struct=fmstruct)
            CALL cp_fm_struct_release(fmstruct=fmstruct)

            CALL cp_fm_to_fm_submat(work%ctransformed(ispin)%matrix, cmatocc, basis_set_iatom%nsgf, nactive(ispin), &
                                    first_sgf(iatom), 1, 1, 1)
            CALL cp_fm_to_fm_submat(ctransformed_virtual(ispin)%matrix, cmatvir, basis_set_iatom%nsgf, nvir(ispin), &
                                    first_sgf(iatom), 1, 1, 1)
!
            CALL cp_gemm(transa="T", &
                         transb="N", &
                         m=nvir(ispin), &
                         n=nactive(ispin), &
                         k=basis_set_iatom%nsgf, &
                         alpha=1.0_dp, &
                         matrix_a=cmatvir, &
                         matrix_b=cmatocc, &
                         beta=0.0_dp, &
                         matrix_c=charges(ispin, iatom)%matrix)

            CALL cp_fm_release(matrix=cmatocc)
            CALL cp_fm_release(matrix=cmatvir)

         END DO !n_atoms
         CALL cp_fm_struct_release(fmstruct=fmstructia)
      END DO !nspins

      NULLIFY (DeltaK)
      ALLOCATE (DeltaK(nspins))

      DO ispin = 1, nspins
         ALLOCATE (charge_per_atom_2index(n_atoms, nactive(ispin)*nvir(ispin)))
         ALLOCATE (gcharge_per_atom_2index(n_atoms, nactive(ispin)*nvir(ispin)))
         charge_per_atom_2index = 0.0_dp
         gcharge_per_atom_2index = 0.0_dp

         IF (debug_this_module) THEN
            ALLOCATE (charge_per_atom(n_atoms))
            charge_per_atom = 0.0_dp
         END IF
         !
         DO iatom = 1, n_atoms
            tmp = 0.0_dp
            DO ivir = 1, nvir(ispin)
               DO iocc = 1, nactive(ispin)
                  kk = (iocc - 1)*nvir(ispin) + ivir
                  CALL cp_fm_get_element(charges(ispin, iatom)%matrix, ivir, iocc, tmp)
                  IF (debug_this_module) charge_per_atom(iatom) = charge_per_atom(iatom) + tmp
                  charge_per_atom_2index(iatom, kk) = charge_per_atom_2index(iatom, kk) + tmp
               END DO
            END DO
         END DO !iatom

         IF ((debug_this_module)) THEN
            ! check charge per atom
            IF (iounit > 0) THEN
             IF ((itertddft .EQ. 1) .AND. (ivect .EQ. 1)) THEN
               DO iatom = 1, n_atoms
                  WRITE (iounit, *) "Charge for Atom", iatom, "   ", charge_per_atom(iatom)
               END DO
             END IF
            END IF
            DEALLOCATE (charge_per_atom)
         END IF

         NULLIFY (fmstructia, tmpmat)
         CALL cp_fm_struct_create(fmstruct=fmstructia, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nvir(ispin), &
                                  ncol_global=nactive(ispin))
         CALL cp_fm_create(matrix=DeltaK(ispin)%matrix, matrix_struct=fmstructia)
         CALL cp_fm_create(matrix=tmpmat, matrix_struct=fmstructia)
         CALL cp_fm_struct_release(fmstruct=fmstructia)

         CALL cp_fm_set_all(DeltaK(ispin)%matrix, 0.0_dp)
         DO iatom = 1, n_atoms
            CALL cp_fm_set_all(gcharges(ispin, iatom)%matrix, 0.0_dp)
         END DO

         tempmat => work%gamma_exchange(1)%matrix
         CALL dbcsr_iterator_start(iter, tempmat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, gab, blk)
            CALL cp_fm_scale_and_add(1.0_dp, gcharges(ispin, iatom)%matrix, gab(1, 1), charges(ispin, jatom)%matrix)
            IF (iatom /= jatom) THEN
               CALL cp_fm_scale_and_add(1.0_dp, gcharges(ispin, jatom)%matrix, gab(1, 1), charges(ispin, iatom)%matrix)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)

         IF (.NOT. ldo_ewald) THEN
            nlim = get_limit(n_atoms, para_env%num_pe, para_env%mepos)
            DO iatom = nlim(1), nlim(2)
               DO jatom = 1, iatom - 1
                  rij = particle_set(iatom)%r - particle_set(jatom)%r
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij(:)**2))
                  IF (dr > 1.e-6_dp) THEN
                     CALL cp_fm_scale_and_add(1.0_dp, gcharges(ispin, iatom)%matrix, 1.0_dp/dr, charges(ispin, jatom)%matrix)
                     CALL cp_fm_scale_and_add(1.0_dp, gcharges(ispin, jatom)%matrix, 1.0_dp/dr, charges(ispin, iatom)%matrix)
                  END IF
               END DO
            END DO
         END IF !ldo_ewald

         gcharge_per_atom_2index = 0.0_dp
         DO iatom = 1, n_atoms
            tmp = 0.0_dp
            DO ivir = 1, nvir(ispin)
               DO iocc = 1, nactive(ispin)
                  kk = (iocc - 1)*nvir(ispin) + ivir
                  CALL cp_fm_get_element(gcharges(ispin, iatom)%matrix, ivir, iocc, tmp)
                  gcharge_per_atom_2index(iatom, kk) = gcharge_per_atom_2index(iatom, kk) + tmp
               END DO
            END DO
         END DO

         IF (ldo_ewald) THEN
            ALLOCATE (charge_per_atom(n_atoms))
            ALLOCATE (gcharge_per_atom(n_atoms, 1))
            CPASSERT(ASSOCIATED(work%ewald_env))
            CPASSERT(ASSOCIATED(work%ewald_pw))
            ewald_env => work%ewald_env
            ewald_pw => work%ewald_pw
            CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
            CALL get_qs_env(qs_env=qs_env, virial=virial, atprop=atprop)
            use_virial = .FALSE.
            calculate_forces = .FALSE.
            n_list => sub_env%sab_orb
            DO iocc = 1, nactive(ispin)
               DO ivir = 1, nvir(ispin)
                  kk = (iocc - 1)*nvir(ispin) + ivir
                  DO iatom = 1, n_atoms
                     gcharge_per_atom(iatom, 1) = gcharge_per_atom_2index(iatom, kk)
                     charge_per_atom(iatom) = charge_per_atom_2index(iatom, kk)
                  END DO
                  CALL tb_ewald_overlap(gcharge_per_atom, charge_per_atom, &
                                        alpha, n_list, virial, use_virial, atprop)
                  CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                                        gcharge_per_atom, charge_per_atom, &
                                        calculate_forces, virial, use_virial, atprop)
                  !self charge interaction contribution
                  IF (para_env%ionode) THEN
                     prefactor = -2._dp*alpha*oorootpi
                     DO iatom = 1, n_atoms
                        gcharge_per_atom_2index(iatom, kk) = gcharge_per_atom(iatom, 1) + prefactor*charge_per_atom(iatom)
                        charge_per_atom_2index(iatom, kk) = charge_per_atom(iatom)
                     END DO
                  END IF
               END DO !ivir
            END DO !iocc
            CALL mp_sum(gcharge_per_atom_2index, para_env%group)
            DEALLOCATE (charge_per_atom, gcharge_per_atom)
         END IF !do_ewald

         DO iocc = 1, nactive(ispin)
            DO ivir = 1, nvir(ispin)
               kk = (iocc - 1)*nvir(ispin) + ivir
               tmp = 0.0_dp
               DO iatom = 1, n_atoms
                  tmp = tmp + gcharge_per_atom_2index(iatom, kk)*charge_per_atom_2index(iatom, kk) ! no spinfac
               END DO
               CALL cp_fm_set_element(DeltaK(ispin)%matrix, ivir, iocc, tmp)
              IF ((itertddft .EQ. 1) .AND. (ivect .EQ. 1)) THEN
               IF (debug_this_module .AND. (iounit > 0)) WRITE (iounit, *) "integral", iocc, ivir, tmp
              ENDIF
            END DO
         END DO

         CALL cp_fm_release(matrix=tmpmat)
         DEALLOCATE (charge_per_atom_2index)
         DEALLOCATE (gcharge_per_atom_2index)
      END DO ! nspins

! calculate X_ia
      DO ispin = 1, nspins
! the shift is evaluated based on (ia|ia)^K
         IF ((shift_delta .GT. 0.0_dp) .AND. (shift_sigma .GT. 0.0_dp)) THEN

            IF ((itertddft .EQ. 1) .AND. (ivect .EQ. 1)) THEN
               IF (lshift_output .AND. iounit > 0) THEN
                  WRITE (iounit, "(/,1X,A,2X,A,6X,A,6X,A)") "Original K shift for transition", "i  ->   a", "a.u.", "eV"
               END IF
            END IF
            DO aa = 1, nvir(ispin)
               DO ii = 1, nactive(ispin)
                  CALL cp_fm_get_element(DeltaK(ispin)%matrix, aa, ii, tmp)
                  steepness = 4.0_dp
                  factor = (tmp/shift_sigma)**steepness + 1.0_dp
                  factortwo = shift_delta/factor
                  CALL cp_fm_set_element(DeltaK(ispin)%matrix, aa, ii, factortwo)
                  IF ((itertddft .EQ. 1) .AND. (ivect .EQ. 1)) THEN
                     IF (lshift_output .AND. iounit > 0) THEN
                        WRITE (iounit, "(30X,I4,2X,A,I4,F10.4,F10.4)") &
                           ii, "->", nactive(ispin) + aa, factortwo, factortwo*evolt
                     END IF
                  END IF
               END DO
            END DO

         ELSE
            CPABORT("Given DIAG_SHIFT values for sTDA not reasonable")
         END IF
      END DO

! Calculate X_ia and then Delta E = \Delta_{ia,ia} X_{ia} ->    via X_ia = X_imu S_munu C_nua
      CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrixkp_s)
      CPASSERT(ASSOCIATED(matrixkp_s))
      IF (SIZE(matrixkp_s, 2) > 1) CPWARN("not implemented for k-points.")
      sm_s => matrixkp_s(1, 1)%matrix

      ! calculate X_ai= C_mua^T S_munu X_nui ; DeltaX_ai = Delta_ai X_ai ; X^{shifted}_mui =C_mua DeltaX_ai
      NULLIFY (SC, Xmo, diagshiftmo, diagshift)
      ALLOCATE (SC(nspins), Xmo(nspins), diagshiftmo(nspins), diagshift(nspins))

      DO ispin = 1, nspins
         NULLIFY (fmstructia, fmstruct)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=SC(ispin)%matrix, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=diagshift(ispin)%matrix, matrix_struct=fmstruct)

         CALL cp_fm_struct_create(fmstruct=fmstructia, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nvir(ispin), &
                                  ncol_global=nactive(ispin))
         CALL cp_fm_create(matrix=Xmo(ispin)%matrix, matrix_struct=fmstructia)
         CALL cp_fm_create(matrix=diagshiftmo(ispin)%matrix, matrix_struct=fmstructia)

         CALL cp_fm_struct_release(fmstruct=fmstructia)

         CALL cp_dbcsr_sm_fm_multiply(sm_s, gs_mos(ispin)%mos_virt, &
                                      SC(ispin)%matrix, nvir(ispin), alpha=1.0_dp, beta=0.0_dp)

         CALL cp_gemm(transa="T", &
                      transb="N", &
                      m=nvir(ispin), &
                      n=nactive(ispin), &
                      k=nao, &
                      alpha=1.0_dp, &
                      matrix_a=SC(ispin)%matrix, &
                      matrix_b=X(ispin)%matrix, &
                      beta=0.0_dp, &
                      matrix_c=Xmo(ispin)%matrix)

         CALL cp_fm_schur_product(DeltaK(ispin)%matrix, Xmo(ispin)%matrix, diagshiftmo(ispin)%matrix)
!
         CALL cp_gemm(transa="N", &
                      transb="N", &
                      m=nao, &
                      n=nactive(ispin), &
                      k=nvir(ispin), &
                      alpha=1.0_dp, &
                      matrix_a=gs_mos(ispin)%mos_virt, &
                      matrix_b=diagshiftmo(ispin)%matrix, &
                      beta=0.0_dp, &
                      matrix_c=diagshift(ispin)%matrix)

         CALL cp_fm_scale_and_add(1.0_dp, res(ispin)%matrix, 1.0_dp, diagshift(ispin)%matrix)

      END DO !nspins

      DO ispin = 1, nspins
         CALL cp_fm_release(matrix=SC(ispin)%matrix)
         CALL cp_fm_release(matrix=diagshift(ispin)%matrix)
         CALL cp_fm_release(matrix=diagshiftmo(ispin)%matrix)
         CALL cp_fm_release(matrix=Xmo(ispin)%matrix)
         CALL cp_fm_release(matrix=DeltaK(ispin)%matrix)
         CALL cp_fm_release(matrix=ctransformed_virtual(ispin)%matrix)
         DO iatom = 1, n_atoms
            CALL cp_fm_release(matrix=charges(ispin, iatom)%matrix)
            CALL cp_fm_release(matrix=gcharges(ispin, iatom)%matrix)
         END DO
      END DO
      DEALLOCATE (SC)
      DEALLOCATE (diagshift)
      DEALLOCATE (diagshiftmo)
      DEALLOCATE (Xmo)
      DEALLOCATE (DeltaK)
      DEALLOCATE (gcharges)
      DEALLOCATE (charges)
      DEALLOCATE (ctransformed_virtual)

      DEALLOCATE (first_sgf)
      DEALLOCATE (last_sgf)
      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE
! **************************************************************************************************
!> \brief ...Calculate the sTDA kernel contribution by contracting the Lowdin MO coefficients --
!>           transition charges with the Coulomb-type or exchange-type integrals
!> \param qs_env ...
!> \param stda_control ...
!> \param stda_env ...
!> \param sub_env ...
!> \param work ...
!> \param is_rks_triplets ...
!> \param X ...
!> \param res ... vector AX with A being the sTDA matrix and X the Davidson trial vector of the
!>                eigenvalue problem A*X = omega*X
!> \param gs_mos ...
!> \param ivect ...
! **************************************************************************************************
   SUBROUTINE stda_calculate_kernel(qs_env, stda_control, stda_env, sub_env, &
                                    work, is_rks_triplets, X, res, gs_mos, ivect)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_control_type)                            :: stda_control
      TYPE(stda_env_type)                                :: stda_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work
      LOGICAL, INTENT(IN)                                :: is_rks_triplets
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: X
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: res
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: ivect

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_calculate_kernel'

      INTEGER                                            :: blk, ewald_type, handle, ia, iatom, &
                                                            ikind, iounit, is, ispin, jatom, &
                                                            jkind, jspin, natom, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, kind_of, last_sgf
      INTEGER, DIMENSION(2)                              :: nactive, nlim
      LOGICAL                                            :: calculate_forces, do_coulomb, do_ewald, &
                                                            do_exchange, ldiagshift, &
                                                            lshift_output, nos_present, use_virial
      REAL(KIND=dp)                                      :: alpha, bp, dr, eta, gabr, hfx, rbeta, &
                                                            spinfac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tcharge, tv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gtcharge
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gab, pblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: xtransformed
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstructjspin
      TYPE(cp_fm_type)                                   :: cvec, cvecjspin
      TYPE(cp_fm_type), POINTER                          :: ct, ctjspin, xt
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: pdens
      TYPE(dbcsr_type), POINTER                          :: tempmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      ldiagshift = stda_control%diag_shift_orig
      lshift_output = stda_control%shift_output

      nactive(:) = stda_env%nactive(:)
      nspins = SIZE(X)
      spinfac = 2.0_dp
      IF (nspins == 2) spinfac = 1.0_dp

      IF (nspins == 1 .AND. is_rks_triplets) THEN
         do_coulomb = .FALSE.
      ELSE
         do_coulomb = .TRUE.
      END IF
      do_ewald = stda_control%do_ewald
      do_exchange = stda_control%do_exchange

      para_env => sub_env%para_env

      CALL get_qs_env(qs_env, natom=natom, cell=cell, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      ! calculate Loewdin transformed Davidson trial vector tilde(X)=S^1/2*X
      ! and tilde(tilde(X))=S^1/2_A*tilde(X)_A
      ALLOCATE (xtransformed(nspins))
      DO ispin = 1, nspins
         NULLIFY (fmstruct)
         ct => work%ctransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct)
         ALLOCATE (xtransformed(ispin)%matrix)
         CALL cp_fm_create(matrix=xtransformed(ispin)%matrix, matrix_struct=fmstruct, name="XTRANSFORMED")
      END DO
      CALL get_lowdin_x(work%shalf, X, xtransformed)

      ALLOCATE (tcharge(natom), gtcharge(natom, 1))
      DO ispin = 1, nspins
         CALL cp_fm_set_all(res(ispin)%matrix, 0.0_dp)
      END DO
      nos_present = .FALSE.

      DO ispin = 1, nspins
         ct => work%ctransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct, nrow_global=nsgf)
         ALLOCATE (tv(nsgf))
         CALL cp_fm_create(cvec, fmstruct)
         !
         ! *** Coulomb contribution
         !
         IF (do_coulomb) THEN
            tcharge(:) = 0.0_dp
            DO jspin = 1, nspins
               ctjspin => work%ctransformed(jspin)%matrix
               xt => xtransformed(jspin)%matrix
               CALL cp_fm_get_info(ctjspin, matrix_struct=fmstructjspin, nrow_global=nsgf)
               CALL cp_fm_create(cvecjspin, fmstructjspin)
               ! CV(mu,j) = CT(mu,j)*XT(mu,j)
               CALL cp_fm_schur_product(ctjspin, xt, cvecjspin)
               ! TV(mu) = SUM_j CV(mu,j)
               CALL cp_fm_vectorssum(cvecjspin, tv, "R")
               ! contract charges
               ! TC(a) = SUM_(mu of a) TV(mu)
               DO ia = 1, natom
                  DO is = first_sgf(ia), last_sgf(ia)
                     tcharge(ia) = tcharge(ia) + tv(is)
                  END DO
               END DO
               CALL cp_fm_release(cvecjspin)
            END DO !jspin

            ! Apply tcharge*gab -> gtcharge
            ! gT(b) = SUM_a g(a,b)*TC(a)
            ! gab = work%gamma_exchange(1)%matrix
            gtcharge = 0.0_dp
            ! short range contribution
            tempmat => work%gamma_exchange(1)%matrix
            CALL dbcsr_iterator_start(iter, tempmat)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, gab, blk)
               gtcharge(iatom, 1) = gtcharge(iatom, 1) + gab(1, 1)*tcharge(jatom)
               IF (iatom /= jatom) THEN
                  gtcharge(jatom, 1) = gtcharge(jatom, 1) + gab(1, 1)*tcharge(iatom)
               END IF
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! Ewald long range contribution
            IF (do_ewald) THEN
               ewald_env => work%ewald_env
               ewald_pw => work%ewald_pw
               CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
               CALL get_qs_env(qs_env=qs_env, virial=virial, atprop=atprop)
               use_virial = .FALSE.
               calculate_forces = .FALSE.
               n_list => sub_env%sab_orb
               CALL tb_ewald_overlap(gtcharge, tcharge, alpha, n_list, virial, use_virial, atprop)
               CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                                     gtcharge, tcharge, calculate_forces, virial, use_virial, atprop)
               ! add self charge interaction contribution
               IF (para_env%ionode) THEN
                  gtcharge(:, 1) = gtcharge(:, 1) - 2._dp*alpha*oorootpi*tcharge(:)
               END IF
            ELSE
               nlim = get_limit(natom, para_env%num_pe, para_env%mepos)
               DO iatom = nlim(1), nlim(2)
                  DO jatom = 1, iatom - 1
                     rij = particle_set(iatom)%r - particle_set(jatom)%r
                     rij = pbc(rij, cell)
                     dr = SQRT(SUM(rij(:)**2))
                     IF (dr > 1.e-6_dp) THEN
                        gtcharge(iatom, 1) = gtcharge(iatom, 1) + tcharge(jatom)/dr
                        gtcharge(jatom, 1) = gtcharge(jatom, 1) + tcharge(iatom)/dr
                     END IF
                  END DO
               END DO
            END IF
            CALL mp_sum(gtcharge, para_env%group)
            ! expand charges
            ! TV(mu) = TC(a of mu)
            tv(1:nsgf) = 0.0_dp
            DO ia = 1, natom
               DO is = first_sgf(ia), last_sgf(ia)
                  tv(is) = gtcharge(ia, 1)
               END DO
            END DO
            ! CV(mu,i) = TV(mu)*CV(mu,i)
            ct => work%ctransformed(ispin)%matrix
            CALL cp_fm_to_fm(ct, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            ! rho(nu,i) = rho(nu,i) + Shalf(nu,mu)*CV(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, res(ispin)%matrix, nactive(ispin), spinfac, 1.0_dp)
         END IF
         !
         ! *** Exchange contribution
         !
         IF (do_exchange) THEN ! option to explicitly switch off exchange
            ! (exchange contributes also if hfx_fraction=0)
            ALLOCATE (kind_of(natom))
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
            !
            tempmat => work%shalf
            CALL dbcsr_create(pdens, template=tempmat, matrix_type=dbcsr_type_no_symmetry)
            ! P(nu,mu) = SUM_j XT(nu,j)*CT(mu,j)
            ct => work%ctransformed(ispin)%matrix
            xt => xtransformed(ispin)%matrix
            CALL dbcsr_set(pdens, 0.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(pdens, xt, ct, nactive(ispin), &
                                       1.0_dp, keep_sparsity=.FALSE.)
            CALL dbcsr_filter(pdens, stda_env%eps_td_filter)
            ! Apply PP*gab -> PP; gab = gamma_coulomb
            ! P(nu,mu) = P(nu,mu)*g(a of nu,b of mu)
            bp = stda_env%beta_param
            hfx = stda_env%hfx_fraction
            CALL dbcsr_iterator_start(iter, pdens)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, pblock, blk)
               rij = particle_set(iatom)%r - particle_set(jatom)%r
               rij = pbc(rij, cell)
               dr = SQRT(SUM(rij(:)**2))
               ikind = kind_of(iatom)
               jkind = kind_of(jatom)
               eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param + &
                      stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp
               rbeta = dr**bp
               IF (hfx > 0.0_dp) THEN
                  gabr = (1._dp/(rbeta + (hfx*eta)**(-bp)))**(1._dp/bp)
               ELSE
                  IF (dr < 1.e-6) THEN
                     gabr = 0.0_dp
                  ELSE
                     gabr = 1._dp/dr
                  END IF
               END IF
               pblock = gabr*pblock
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! CV(mu,i) = P(nu,mu)*CT(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(pdens, ct, cvec, nactive(ispin), 1.0_dp, 0.0_dp)
            ! rho(nu,i) = rho(nu,i) + ShalfP(nu,mu)*CV(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, res(ispin)%matrix, nactive(ispin), -1.0_dp, 1.0_dp)
            !
            CALL dbcsr_release(pdens)
            DEALLOCATE (kind_of)
         END IF
         !
         CALL cp_fm_release(cvec)
         DEALLOCATE (tv)
      END DO

      IF (ldiagshift) THEN
         CALL k_shift_orig(qs_env, sub_env, stda_control, work, gs_mos, X, ivect, res)
      END IF
      IF (lshift_output) THEN
         CALL more_output(qs_env, sub_env, stda_control, work, gs_mos, ivect)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_release(matrix=xtransformed(ispin)%matrix)
         DEALLOCATE (xtransformed(ispin)%matrix)
      END DO
      DEALLOCATE (xtransformed)
      DEALLOCATE (tcharge, gtcharge)
      DEALLOCATE (first_sgf, last_sgf)

      CALL timestop(handle)

   END SUBROUTINE stda_calculate_kernel

! **************************************************************************************************
END MODULE qs_tddfpt2_stda_utils
