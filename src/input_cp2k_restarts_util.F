!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE input_cp2k_restarts_util
   USE cp_linked_list_input,            ONLY: cp_sll_val_create,&
                                              cp_sll_val_get_length,&
                                              cp_sll_val_type
   USE input_section_types,             ONLY: section_get_keyword_index,&
                                              section_type,&
                                              section_vals_add_values,&
                                              section_vals_type
   USE input_val_types,                 ONLY: val_create,&
                                              val_release,&
                                              val_type
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_restarts_util'

   PUBLIC :: section_velocity_val_set

CONTAINS

! **************************************************************************************************
!> \brief routine to dump velocities.. fast implementation
!> \param velocity_section ...
!> \param particles ...
!> \param velocity ...
!> \param conv_factor ...
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE section_velocity_val_set(velocity_section, particles, velocity, conv_factor)

      TYPE(section_vals_type), POINTER                   :: velocity_section
      TYPE(particle_list_type), OPTIONAL, POINTER        :: particles
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: velocity
      REAL(KIND=dp)                                      :: conv_factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'section_velocity_val_set'

      INTEGER                                            :: handle, ik, irk, Nlist, nloop
      LOGICAL                                            :: check
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vel
      TYPE(cp_sll_val_type), POINTER                     :: new_pos, vals
      TYPE(section_type), POINTER                        :: section
      TYPE(val_type), POINTER                            :: my_val, old_val

      CALL timeset(routineN, handle)
      NULLIFY (my_val, old_val, section, vals)
      CPASSERT(ASSOCIATED(velocity_section))
      CPASSERT(velocity_section%ref_count > 0)
      section => velocity_section%section
      ik = section_get_keyword_index(section, "_DEFAULT_KEYWORD_")
      IF (ik == -2) &
         CALL cp_abort(__LOCATION__, &
                       "section "//TRIM(section%name)//" does not contain keyword "// &
                       "_DEFAULT_KEYWORD_")

      ! At least one of the two arguments must be present..
      check = PRESENT(particles) .NEQV. PRESENT(velocity)
      CPASSERT(check)
      IF (PRESENT(particles)) nloop = particles%n_els
      IF (PRESENT(velocity)) nloop = SIZE(velocity, 2)

      DO
         IF (SIZE(velocity_section%values, 2) == 1) EXIT
         CALL section_vals_add_values(velocity_section)
      END DO
      vals => velocity_section%values(ik, 1)%list
      Nlist = 0
      IF (ASSOCIATED(vals)) THEN
         Nlist = cp_sll_val_get_length(vals)
      END IF
      DO irk = 1, nloop
         ALLOCATE (vel(3))
         ! Always stored in A.U.
         IF (PRESENT(particles)) vel = particles%els(irk)%v(1:3)*conv_factor
         IF (PRESENT(velocity)) vel = velocity(1:3, irk)*conv_factor

         CALL val_create(my_val, r_vals_ptr=vel)

         IF (Nlist /= 0) THEN
            IF (irk == 1) THEN
               new_pos => vals
            ELSE
               new_pos => new_pos%rest
            END IF
            old_val => new_pos%first_el
            CALL val_release(old_val)
            new_pos%first_el => my_val
         ELSE
            IF (irk == 1) THEN
               NULLIFY (new_pos)
               CALL cp_sll_val_create(new_pos, first_el=my_val)
               vals => new_pos
            ELSE
               NULLIFY (new_pos%rest)
               CALL cp_sll_val_create(new_pos%rest, first_el=my_val)
               new_pos => new_pos%rest
            END IF
         END IF
         NULLIFY (my_val)
      END DO

      velocity_section%values(ik, 1)%list => vals

      CALL timestop(handle)

   END SUBROUTINE section_velocity_val_set

END MODULE input_cp2k_restarts_util
