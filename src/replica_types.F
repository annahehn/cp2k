!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types used to handle many replica of the same system that differ only
!>      in atom positions, and velocity.
!>      This is useful for things like path integrals or nudged elastic band
!> \note
!>      this is a stupid implementation that replicates all the information
!>      about the replicas, if you really want to do a *lot* of replicas on
!>      a lot of processors you should think about distributiong also that
!>      information
!> \par History
!>      09.2005 created [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE replica_types
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_rm_iter_level
   USE cp_result_methods,               ONLY: cp_results_mp_bcast
   USE cp_result_types,                 ONLY: cp_result_p_type,&
                                              cp_result_release
   USE f77_interface,                   ONLY: destroy_force_env,&
                                              f_env_add_defaults,&
                                              f_env_rm_defaults,&
                                              f_env_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_cart_release,&
                                              mp_para_cart_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE qs_wf_history_types,             ONLY: qs_wf_history_p_type,&
                                              wfi_release
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   LOGICAL, SAVE, PRIVATE :: module_initialized = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_types'

   PUBLIC :: replica_env_type
   PUBLIC :: rep_env_release
   PUBLIC :: rep_env_sync, rep_env_sync_results, rep_envs_add_rep_env
   PUBLIC :: rep_envs_get_rep_env

! **************************************************************************************************
!> \brief keeps replicated information about the replicas
!> \param ref_count reference count
!> \param id_nr identity number (unique or each replica_env)
!> \param nrep number of replicas
!> \param nparticle number of particles (usually atoms) in each replica
!> \param ndim = 3*nparticle
!> \param f_env_id id of the force env that will do the calculations for the
!>        replicas owned by this processor
!> \param r ,v,f: positions, velocities and forces of the replicas.
!>        the indexing is as follow (idir,iat,irep)
!> \param replica_owner which replica group number owns the replica irep
!> \param cart 2d distribution of the processors for the replicas,
!>        a column (or row if row_force was true in the rep_env_create call)
!>        work together on the same force_env (i.e. changing the
!>        row (column) you stay in the same replica), rows (columns) have
!>        different replicas
!> \param force_dim which dimension of cart works on forces together
!>        used to be hardcoded to 1. Default is still 1, will
!>        be 2 if row_force is true in the rep_env_create call.
!> \param para_env the global para env that contains all the replicas,
!>        this is just the cart as para_env
!> \param para_env_f parallel environment of the underlying force
!>        environment
!> \param inter_rep_rank mapping replica group number -> rank in para_env_inter_rep
!>        (this used to be col_rank)
!> \param para_env_inter_rep parallel environment between replica
!> \param force_rank mapping number of processor in force env -> rank in para_env_f
!>        (this used to be row_rank)
!> \param local_rep_indices indices of the local replicas, starting at 1
!> \param rep_is_local logical if specific replica is a local one.
!> \param my_rep_group which replica group number this process belongs to
!>        (this used to be just cart%mepos(2) but with transposing the cart
!>        (row_force=.true.) became cart%mepos(1), and to generalize this it
!>        is now a separate variable, so one does not need to know
!>        which way the cart is mapped.)
!> \param wf_history wavefunction history for the owned replicas
!> \param keep_wf_history if the wavefunction history for the owned replicas
!>        should be kept
!> \author fawzi
! **************************************************************************************************
   TYPE replica_env_type
      INTEGER                                           :: ref_count = -1, id_nr = -1, f_env_id = -1, &
                                                           nrep = -1, ndim = -1, nparticle = -1, &
                                                           my_rep_group = -1, force_dim = -1
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: r => NULL(), v => NULL(), f => NULL()
      LOGICAL                                           :: sync_v = .FALSE., keep_wf_history = .FALSE.
      CHARACTER(LEN=default_path_length)                :: original_project_name = ""
      TYPE(qs_wf_history_p_type), DIMENSION(:), POINTER :: wf_history => NULL()
      TYPE(cp_result_p_type), DIMENSION(:), POINTER       :: results => NULL()
      INTEGER, DIMENSION(:), POINTER                    :: local_rep_indices => NULL()
      INTEGER, DIMENSION(:), POINTER                    :: replica_owner => NULL(), force_rank => NULL(), &
                                                           inter_rep_rank => NULL()
      LOGICAL, DIMENSION(:), POINTER                    :: rep_is_local => NULL()
      TYPE(mp_para_cart_type), POINTER                  :: cart => NULL()
      TYPE(mp_para_env_type), POINTER                   :: para_env => NULL(), para_env_f => NULL(), &
                                                           para_env_inter_rep => NULL()
   END TYPE replica_env_type

! **************************************************************************************************
!> \brief ****s* replica_types/replica_env_p_type *
!>
!>      to build arrays of pointers to a replica_env_type
!> \param rep_env the pointer to the replica_env
!> \author fawzi
! **************************************************************************************************
   TYPE replica_env_p_type
      TYPE(replica_env_type), POINTER                   :: rep_env => NULL()
   END TYPE replica_env_p_type

   TYPE(replica_env_p_type), POINTER, DIMENSION(:), PRIVATE :: rep_envs

CONTAINS

! **************************************************************************************************
!> \brief releases the given replica environment
!> \param rep_env the replica environment to release
!> \author fawzi
!> \note
!>      here and not in replica_types to allow the use of replica_env_type
!>      in a force_env (call to destroy_force_env gives circular dep)
! **************************************************************************************************
   SUBROUTINE rep_env_release(rep_env)
      TYPE(replica_env_type), POINTER                    :: rep_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rep_env_release'

      INTEGER                                            :: handle, i, ierr

      CALL timeset(routineN, handle)
      IF (ASSOCIATED(rep_env)) THEN
         CPASSERT(rep_env%ref_count > 0)
         rep_env%ref_count = rep_env%ref_count - 1
         IF (rep_env%ref_count == 0) THEN
            CALL rep_env_destroy_low(rep_env%id_nr, ierr)
            IF (rep_env%f_env_id > 0) THEN
               CALL destroy_force_env(rep_env%f_env_id, ierr)
               CPASSERT(ierr == 0)
            END IF
            IF (ASSOCIATED(rep_env%r)) THEN
               DEALLOCATE (rep_env%r)
            END IF
            IF (ASSOCIATED(rep_env%v)) THEN
               DEALLOCATE (rep_env%v)
            END IF
            IF (ASSOCIATED(rep_env%f)) THEN
               DEALLOCATE (rep_env%f)
            END IF
            IF (ASSOCIATED(rep_env%wf_history)) THEN
               DO i = 1, SIZE(rep_env%wf_history)
                  CALL wfi_release(rep_env%wf_history(i)%wf_history)
               END DO
               DEALLOCATE (rep_env%wf_history)
            END IF
            IF (ASSOCIATED(rep_env%results)) THEN
               DO i = 1, SIZE(rep_env%results)
                  CALL cp_result_release(rep_env%results(i)%results)
               END DO
               DEALLOCATE (rep_env%results)
            END IF
            DEALLOCATE (rep_env%local_rep_indices)
            DEALLOCATE (rep_env%rep_is_local)
            IF (ASSOCIATED(rep_env%replica_owner)) THEN
               DEALLOCATE (rep_env%replica_owner)
            END IF
            DEALLOCATE (rep_env%inter_rep_rank, rep_env%force_rank)
            CALL mp_para_cart_release(rep_env%cart)
            CALL mp_para_env_release(rep_env%para_env)
            CALL mp_para_env_release(rep_env%para_env_f)
            CALL mp_para_env_release(rep_env%para_env_inter_rep)
            CALL rep_envs_rm_rep_env(rep_env)
            DEALLOCATE (rep_env)
         END IF
      END IF
      NULLIFY (rep_env)
      CALL timestop(handle)
   END SUBROUTINE rep_env_release

! **************************************************************************************************
!> \brief initializes the destruction of the replica_env
!> \param rep_env_id id_nr of the replica environment that should be initialized
!> \param ierr will be non zero if there is an initialization error
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE rep_env_destroy_low(rep_env_id, ierr)
      INTEGER, INTENT(in)                                :: rep_env_id
      INTEGER, INTENT(out)                               :: ierr

      INTEGER                                            :: stat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(replica_env_type), POINTER                    :: rep_env

      rep_env => rep_envs_get_rep_env(rep_env_id, ierr=stat)
      IF (.NOT. ASSOCIATED(rep_env)) &
         CPABORT("could not find rep_env with id_nr"//cp_to_string(rep_env_id))
      CALL f_env_add_defaults(f_env_id=rep_env%f_env_id, f_env=f_env)
      logger => cp_get_default_logger()
      CALL cp_rm_iter_level(iteration_info=logger%iter_info, &
                            level_name="REPLICA_EVAL")
      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE rep_env_destroy_low

! **************************************************************************************************
!> \brief sends the data from each replica to all the other
!>      on replica j/=i data from replica i overwrites val(:,i)
!> \param rep_env replica environment
!> \param vals the values to synchronize (second index runs over replicas)
!> \author fawzi
!> \note
!>      could be optimized: bcast in inter_rep, all2all or shift vs sum
! **************************************************************************************************
   SUBROUTINE rep_env_sync(rep_env, vals)
      TYPE(replica_env_type), POINTER                    :: rep_env
      REAL(kind=dp), DIMENSION(:, :), INTENT(inout)      :: vals

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rep_env_sync'

      INTEGER                                            :: handle, irep

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(rep_env))
      CPASSERT(rep_env%ref_count > 0)
      CPASSERT(SIZE(vals, 2) == rep_env%nrep)
      DO irep = 1, rep_env%nrep
         IF (.NOT. rep_env%rep_is_local(irep)) THEN
            vals(:, irep) = 0._dp
         END IF
      END DO
      CALL rep_env%para_env_inter_rep%sum(vals)
      CALL timestop(handle)
   END SUBROUTINE rep_env_sync

! **************************************************************************************************
!> \brief sends the data from each replica to all the other
!>      in this case the result type is passed
!> \param rep_env replica environment
!> \param results is an array of result_types
!> \author fschiff
! **************************************************************************************************
   SUBROUTINE rep_env_sync_results(rep_env, results)
      TYPE(replica_env_type), POINTER                    :: rep_env
      TYPE(cp_result_p_type), DIMENSION(:), POINTER      :: results

      CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_sync_results'

      INTEGER                                            :: handle, irep, nrep, source

      CALL timeset(routineN, handle)
      nrep = rep_env%nrep
      CPASSERT(ASSOCIATED(rep_env))
      CPASSERT(rep_env%ref_count > 0)
      CPASSERT(SIZE(results) == rep_env%nrep)
      DO irep = 1, nrep
         source = rep_env%inter_rep_rank(rep_env%replica_owner(irep))
         CALL cp_results_mp_bcast(results(irep)%results, source, rep_env%para_env_inter_rep)
      END DO
      CALL timestop(handle)
   END SUBROUTINE rep_env_sync_results

! **************************************************************************************************
!> \brief returns the replica environment with the given id_nr
!> \param id_nr the id_nr of the requested rep_envs
!> \param ierr ...
!> \return ...
!> \author fawzi
! **************************************************************************************************
   FUNCTION rep_envs_get_rep_env(id_nr, ierr) RESULT(res)
      INTEGER, INTENT(in)                                :: id_nr
      INTEGER, INTENT(OUT)                               :: ierr
      TYPE(replica_env_type), POINTER                    :: res

      INTEGER                                            :: i

      NULLIFY (res)
      ierr = -1
      IF (module_initialized) THEN
         IF (ASSOCIATED(rep_envs)) THEN
            DO i = 1, SIZE(rep_envs)
               IF (rep_envs(i)%rep_env%id_nr == id_nr) THEN
                  res => rep_envs(i)%rep_env
                  ierr = 0
                  EXIT
               END IF
            END DO
         END IF
      END IF
   END FUNCTION rep_envs_get_rep_env

! **************************************************************************************************
!> \brief adds the given rep_env to the list of controlled rep_envs.
!> \param rep_env the rep_env to add
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE rep_envs_add_rep_env(rep_env)
      TYPE(replica_env_type), POINTER                    :: rep_env

      INTEGER                                            :: i, stat
      TYPE(replica_env_p_type), DIMENSION(:), POINTER    :: new_rep_envs
      TYPE(replica_env_type), POINTER                    :: rep_env2

      IF (ASSOCIATED(rep_env)) THEN
         rep_env2 => rep_envs_get_rep_env(rep_env%id_nr, ierr=stat)
         IF (.NOT. ASSOCIATED(rep_env2)) THEN
            IF (module_initialized) THEN
               IF (.NOT. ASSOCIATED(rep_envs)) THEN
                  ALLOCATE (rep_envs(1))
               ELSE
                  ALLOCATE (new_rep_envs(SIZE(rep_envs) + 1))
                  DO i = 1, SIZE(rep_envs)
                     new_rep_envs(i)%rep_env => rep_envs(i)%rep_env
                  END DO
                  DEALLOCATE (rep_envs)
                  rep_envs => new_rep_envs
               END IF
            ELSE
               ALLOCATE (rep_envs(1))
            END IF
            rep_envs(SIZE(rep_envs))%rep_env => rep_env
            module_initialized = .TRUE.
         END IF
      END IF
   END SUBROUTINE rep_envs_add_rep_env

! **************************************************************************************************
!> \brief removes the given rep_env to the list of controlled rep_envs.
!> \param rep_env the rep_env to remove
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE rep_envs_rm_rep_env(rep_env)
      TYPE(replica_env_type), POINTER                    :: rep_env

      INTEGER                                            :: i, ii
      TYPE(replica_env_p_type), DIMENSION(:), POINTER    :: new_rep_envs

      IF (ASSOCIATED(rep_env)) THEN
         CPASSERT(module_initialized)
         ALLOCATE (new_rep_envs(SIZE(rep_envs) - 1))
         ii = 0
         DO i = 1, SIZE(rep_envs)
            IF (rep_envs(i)%rep_env%id_nr /= rep_env%id_nr) THEN
               ii = ii + 1
               new_rep_envs(ii)%rep_env => rep_envs(i)%rep_env
            END IF
         END DO
         CPASSERT(ii == SIZE(new_rep_envs))
         DEALLOCATE (rep_envs)
         rep_envs => new_rep_envs
         IF (SIZE(rep_envs) == 0) THEN
            DEALLOCATE (rep_envs)
         END IF
      END IF
   END SUBROUTINE rep_envs_rm_rep_env

END MODULE replica_types
