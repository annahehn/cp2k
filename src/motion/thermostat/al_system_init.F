!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author Noam Bernstein [noamb] 02.2012
! **************************************************************************************************
MODULE al_system_init

   USE al_system_mapping,               ONLY: al_to_particle_mapping
   USE al_system_types,                 ONLY: al_system_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: global_constraint_type,&
                                              molecule_type
   USE simpar_types,                    ONLY: simpar_type
   USE thermostat_types,                ONLY: thermostat_info_type
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: initialize_al_part

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'al_system_init'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule ...
!> \param molecule_kind_set ...
!> \param para_env ...
!> \param al ...
!> \param al_section ...
!> \param gci ...
!> \author Noam Bernstein [noamb] 02.2012
! **************************************************************************************************
   SUBROUTINE initialize_al_part(thermostat_info, simpar, local_molecules, &
                                 molecule, molecule_kind_set, para_env, al, al_section, &
                                 gci)

      TYPE(thermostat_info_type), POINTER                :: thermostat_info
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_type), POINTER                       :: molecule(:)
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(al_system_type), POINTER                      :: al
      TYPE(section_vals_type), POINTER                   :: al_section
      TYPE(global_constraint_type), POINTER              :: gci

      LOGICAL                                            :: restart

      restart = .FALSE.
      CALL al_to_particle_mapping(thermostat_info, simpar, local_molecules, &
                                  molecule, molecule_kind_set, al, para_env, gci)

      CALL restart_al(al, al_section, restart)

      IF (.NOT. restart) THEN
         CALL init_al_variables(al)
      END IF

   END SUBROUTINE initialize_al_part

! **************************************************************************************************
!> \brief ...
!> \param al ...
! **************************************************************************************************
   SUBROUTINE init_al_variables(al)
      TYPE(al_system_type), POINTER                      :: al

      al%nvt(:)%mass = al%nvt(:)%nkt*al%tau_nh**2
      al%nvt(:)%chi = 0.0_dp
   END SUBROUTINE init_al_variables

! **************************************************************************************************
!> \brief ...
!> \param al ...
!> \param al_section ...
!> \param restart ...
!> \author Noam Bernstein [noamb] 02.2012
! **************************************************************************************************
   SUBROUTINE restart_al(al, al_section, restart)
      TYPE(al_system_type), POINTER                      :: al
      TYPE(section_vals_type), POINTER                   :: al_section
      LOGICAL, INTENT(inout)                             :: restart

      INTEGER                                            :: i, my_index, n_rep
      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: work_section

      restart = .FALSE.

      ! Possibly restart the initial thermostat DOF value
      work_section => section_vals_get_subs_vals(section_vals=al_section, &
                                                 subsection_name="CHI")
      CALL section_vals_get(work_section, explicit=explicit)
      restart = explicit
      IF (explicit) THEN
         CALL section_vals_val_get(section_vals=work_section, keyword_name="_DEFAULT_KEYWORD_", &
                                   n_rep_val=n_rep)
         IF (n_rep == al%glob_num_al) THEN
            DO i = 1, al%loc_num_al
               my_index = al%map_info%index(i)
               CALL section_vals_val_get(section_vals=work_section, keyword_name="_DEFAULT_KEYWORD_", &
                                         i_rep_val=my_index, r_val=al%nvt(i)%chi)
            END DO
         ELSE
            CALL cp_abort(__LOCATION__, &
                          'Number pf restartable stream not equal to the number of'// &
                          ' total thermostats!')
         END IF
      END IF

      ! Possibly restart the initial thermostat mass
      work_section => section_vals_get_subs_vals(section_vals=al_section, &
                                                 subsection_name="MASS")
      CALL section_vals_get(work_section, explicit=explicit)
      IF (restart .NEQV. explicit) &
         CALL cp_abort(__LOCATION__, &
                       "You need to define both CHI and MASS sections (or none) in the AD_LANGEVIN section")
      restart = restart .AND. explicit
      IF (explicit) THEN
         CALL section_vals_val_get(section_vals=work_section, keyword_name="_DEFAULT_KEYWORD_", &
                                   n_rep_val=n_rep)
         IF (n_rep == al%glob_num_al) THEN
            DO i = 1, al%loc_num_al
               my_index = al%map_info%index(i)
               CALL section_vals_val_get(section_vals=work_section, keyword_name="_DEFAULT_KEYWORD_", &
                                         i_rep_val=my_index, r_val=al%nvt(i)%mass)
            END DO
         ELSE
            CALL cp_abort(__LOCATION__, &
                          'Number pf restartable stream not equal to the number of'// &
                          ' total thermostats!')
         END IF
      END IF

   END SUBROUTINE restart_al

END MODULE al_system_init
