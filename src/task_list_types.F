!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types for task lists
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE task_list_types
   USE grid_api,                        ONLY: grid_basis_set_type,&
                                              grid_free_basis_set,&
                                              grid_free_task_list,&
                                              grid_task_list_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE offload_api,                     ONLY: offload_buffer_type,&
                                              offload_free_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_types'

   TYPE task_type
      INTEGER                                       :: destination = 0
      INTEGER                                       :: source = 0
      INTEGER                                       :: dist_type = 0
      INTEGER                                       :: cost = 0
      INTEGER                                       :: subpatch_pattern = 0
      INTEGER                                       :: grid_level = 0
      INTEGER                                       :: image = 0
      INTEGER                                       :: iatom = 0
      INTEGER                                       :: jatom = 0
      INTEGER                                       :: iset = 0
      INTEGER                                       :: jset = 0
      INTEGER                                       :: ipgf = 0
      INTEGER                                       :: jpgf = 0
      INTEGER                                       :: pair_index = 0

      REAL(KIND=dp), DIMENSION(3)                   :: rab = 0.0_dp
      REAL(KIND=dp)                                 :: radius = 0.0_dp
   END TYPE task_type

   TYPE atom_pair_type
      INTEGER                                       :: rank = 0
      INTEGER                                       :: row = 0
      INTEGER                                       :: col = 0
      INTEGER                                       :: image = 0
   END TYPE atom_pair_type

   TYPE task_list_type
      TYPE(task_type), DIMENSION(:), POINTER        :: tasks => Null()
      INTEGER                                       :: ntasks = 0
      INTEGER, DIMENSION(:, :), POINTER             :: taskstart => Null(), taskstop => Null()
      INTEGER, DIMENSION(:), POINTER                :: npairs => Null()

      TYPE(atom_pair_type), DIMENSION(:), POINTER   :: atom_pair_send => Null()
      INTEGER, DIMENSION(:), POINTER                :: pair_offsets_send => Null()
      INTEGER, DIMENSION(:), POINTER                :: rank_offsets_send => Null()
      INTEGER, DIMENSION(:), POINTER                :: rank_sizes_send => Null()
      INTEGER                                       :: buffer_size_send = 0

      TYPE(atom_pair_type), DIMENSION(:), POINTER   :: atom_pair_recv => Null()
      INTEGER, DIMENSION(:), POINTER                :: pair_offsets_recv => Null()
      INTEGER, DIMENSION(:), POINTER                :: rank_offsets_recv => Null()
      INTEGER, DIMENSION(:), POINTER                :: rank_sizes_recv => Null()
      INTEGER                                       :: buffer_size_recv = 0

      TYPE(grid_basis_set_type), DIMENSION(:), POINTER :: grid_basis_sets => Null()
      TYPE(grid_task_list_type)                     :: grid_task_list = grid_task_list_type()
      TYPE(offload_buffer_type)                     :: pab_buffer = offload_buffer_type()
      TYPE(offload_buffer_type)                     :: hab_buffer = offload_buffer_type()
   END TYPE task_list_type

   INTEGER, PARAMETER                               :: task_size_in_int8 = 17

   PUBLIC :: task_type, atom_pair_type, task_list_type
   PUBLIC :: task_size_in_int8, serialize_task, deserialize_task
   PUBLIC :: allocate_task_list, deallocate_task_list, reallocate_tasks

CONTAINS

! **************************************************************************************************
!> \brief allocates and initialised the components of the task_list_type
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE allocate_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      ALLOCATE (task_list)

      NULLIFY (task_list%tasks)
      NULLIFY (task_list%atom_pair_send)
      NULLIFY (task_list%atom_pair_recv)
      NULLIFY (task_list%taskstart)
      NULLIFY (task_list%taskstop)
      NULLIFY (task_list%npairs)
      task_list%ntasks = 0
   END SUBROUTINE allocate_task_list

! **************************************************************************************************
!> \brief Grow an array of tasks while preserving the existing entries.
!> \param tasks ...
!> \param new_size ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE reallocate_tasks(tasks, new_size)
      TYPE(task_type), DIMENSION(:), INTENT(inout), &
         POINTER                                         :: tasks
      INTEGER, INTENT(in)                                :: new_size

      TYPE(task_type), DIMENSION(:), POINTER             :: old_tasks

      IF (.NOT. ASSOCIATED(tasks)) THEN
         ALLOCATE (tasks(new_size))
         RETURN
      END IF

      IF (new_size <= SIZE(tasks)) RETURN

      old_tasks => tasks
      NULLIFY (tasks)

      ALLOCATE (tasks(new_size))
      tasks(1:SIZE(old_tasks)) = old_tasks(:)
      DEALLOCATE (old_tasks)
   END SUBROUTINE reallocate_tasks

! **************************************************************************************************
!> \brief deallocates the components and the object itself
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE deallocate_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      INTEGER                                            :: i

      IF (ASSOCIATED(task_list%tasks)) THEN
         DEALLOCATE (task_list%tasks)
      END IF
      IF (ASSOCIATED(task_list%taskstart)) THEN
         DEALLOCATE (task_list%taskstart)
      END IF
      IF (ASSOCIATED(task_list%taskstop)) THEN
         DEALLOCATE (task_list%taskstop)
      END IF
      IF (ASSOCIATED(task_list%npairs)) THEN
         DEALLOCATE (task_list%npairs)
      END IF

      IF (ASSOCIATED(task_list%atom_pair_send)) THEN
         DEALLOCATE (task_list%atom_pair_send)
      END IF
      IF (ASSOCIATED(task_list%pair_offsets_send)) THEN
         DEALLOCATE (task_list%pair_offsets_send)
      END IF
      IF (ASSOCIATED(task_list%rank_offsets_send)) THEN
         DEALLOCATE (task_list%rank_offsets_send)
      END IF
      IF (ASSOCIATED(task_list%rank_sizes_send)) THEN
         DEALLOCATE (task_list%rank_sizes_send)
      END IF

      IF (ASSOCIATED(task_list%atom_pair_recv)) THEN
         DEALLOCATE (task_list%atom_pair_recv)
      END IF
      IF (ASSOCIATED(task_list%pair_offsets_recv)) THEN
         DEALLOCATE (task_list%pair_offsets_recv)
      END IF
      IF (ASSOCIATED(task_list%rank_offsets_recv)) THEN
         DEALLOCATE (task_list%rank_offsets_recv)
      END IF
      IF (ASSOCIATED(task_list%rank_sizes_recv)) THEN
         DEALLOCATE (task_list%rank_sizes_recv)
      END IF
      CALL grid_free_task_list(task_list%grid_task_list)
      CALL offload_free_buffer(task_list%pab_buffer)
      CALL offload_free_buffer(task_list%hab_buffer)
      IF (ASSOCIATED(task_list%grid_basis_sets)) THEN
         DO i = 1, SIZE(task_list%grid_basis_sets)
            CALL grid_free_basis_set(task_list%grid_basis_sets(i))
         END DO
         DEALLOCATE (task_list%grid_basis_sets)
      END IF

      DEALLOCATE (task_list)
   END SUBROUTINE deallocate_task_list

! **************************************************************************************************
!> \brief Serialize a task into an integer array. Used for MPI communication.
!> \param task ...
!> \param serialized_task ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE serialize_task(task, serialized_task)
      TYPE(task_type), INTENT(IN)                        :: task
      INTEGER(KIND=int_8), &
         DIMENSION(task_size_in_int8), INTENT(OUT)       :: serialized_task

      serialized_task(1) = task%destination
      serialized_task(2) = task%source
      serialized_task(3) = task%dist_type
      serialized_task(4) = task%cost
      serialized_task(5) = task%subpatch_pattern

      serialized_task(6) = task%grid_level
      serialized_task(7) = task%image
      serialized_task(8) = task%iatom
      serialized_task(9) = task%jatom
      serialized_task(10) = task%iset
      serialized_task(11) = task%jset
      serialized_task(12) = task%ipgf
      serialized_task(13) = task%jpgf

      serialized_task(14) = TRANSFER(task%rab(1), mold=1_int_8)
      serialized_task(15) = TRANSFER(task%rab(2), mold=1_int_8)
      serialized_task(16) = TRANSFER(task%rab(3), mold=1_int_8)
      serialized_task(17) = TRANSFER(task%radius, mold=1_int_8)
   END SUBROUTINE serialize_task

! **************************************************************************************************
!> \brief De-serialize a task from an integer array. Used for MPI communication.
!> \param task ...
!> \param serialized_task ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE deserialize_task(task, serialized_task)
      TYPE(task_type), INTENT(OUT)                       :: task
      INTEGER(KIND=int_8), &
         DIMENSION(task_size_in_int8), INTENT(IN)        :: serialized_task

      task%destination = INT(serialized_task(1))
      task%source = INT(serialized_task(2))
      task%dist_type = INT(serialized_task(3))
      task%cost = INT(serialized_task(4))
      task%subpatch_pattern = INT(serialized_task(5))

      task%grid_level = INT(serialized_task(6))
      task%image = INT(serialized_task(7))
      task%iatom = INT(serialized_task(8))
      task%jatom = INT(serialized_task(9))
      task%iset = INT(serialized_task(10))
      task%jset = INT(serialized_task(11))
      task%ipgf = INT(serialized_task(12))
      task%jpgf = INT(serialized_task(13))

      task%rab(1) = TRANSFER(serialized_task(14), mold=1.0_dp)
      task%rab(2) = TRANSFER(serialized_task(15), mold=1.0_dp)
      task%rab(3) = TRANSFER(serialized_task(16), mold=1.0_dp)
      task%radius = TRANSFER(serialized_task(17), mold=1.0_dp)
   END SUBROUTINE deserialize_task

END MODULE task_list_types
