!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for Thermostats
!> \author teo [tlaino] - University of Zurich - 10.2007
! **************************************************************************************************
MODULE thermostat_methods
   USE al_system_dynamics,              ONLY: al_particles
   USE al_system_init,                  ONLY: initialize_al_part
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE bibliography,                    ONLY: VandeVondele2002,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE csvr_system_dynamics,            ONLY: csvr_barostat,&
                                              csvr_particles,&
                                              csvr_shells
   USE csvr_system_init,                ONLY: initialize_csvr_baro,&
                                              initialize_csvr_part,&
                                              initialize_csvr_shell
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE extended_system_dynamics,        ONLY: lnhc_barostat,&
                                              lnhc_particles,&
                                              lnhc_shells
   USE extended_system_init,            ONLY: initialize_nhc_baro,&
                                              initialize_nhc_fast,&
                                              initialize_nhc_part,&
                                              initialize_nhc_shell,&
                                              initialize_nhc_slow
   USE extended_system_types,           ONLY: npt_info_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE gle_system_dynamics,             ONLY: gle_particles,&
                                              initialize_gle_part
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: &
        do_region_global, do_thermo_al, do_thermo_csvr, do_thermo_gle, do_thermo_nose, &
        do_thermo_same_as_part, npe_f_ensemble, npe_i_ensemble, npt_f_ensemble, npt_i_ensemble, &
        npt_ia_ensemble, nve_ensemble, nvt_adiabatic_ensemble, nvt_ensemble
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_remove_values,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: global_constraint_type,&
                                              molecule_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE qmmm_types,                      ONLY: qmmm_env_type
   USE simpar_types,                    ONLY: simpar_type
   USE thermostat_types,                ONLY: allocate_thermostats,&
                                              create_thermostat_type,&
                                              release_thermostat_info,&
                                              release_thermostat_type,&
                                              release_thermostats,&
                                              thermostat_type,&
                                              thermostats_type
   USE thermostat_utils,                ONLY: compute_degrees_of_freedom,&
                                              compute_nfree,&
                                              get_thermostat_energies,&
                                              setup_adiabatic_thermostat_info,&
                                              setup_thermostat_info
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: create_thermostats, &
             apply_thermostat_baro, &
             apply_thermostat_particles, &
             apply_thermostat_shells

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_methods'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param thermostats ...
!> \param md_section ...
!> \param force_env ...
!> \param simpar ...
!> \param para_env ...
!> \param globenv ...
!> \param global_section ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE create_thermostats(thermostats, md_section, force_env, simpar, &
                                 para_env, globenv, global_section)
      TYPE(thermostats_type), POINTER                    :: thermostats
      TYPE(section_vals_type), POINTER                   :: md_section
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: global_section

      CHARACTER(LEN=default_path_length)                 :: binary_restart_file_name
      INTEGER                                            :: n_rep, region, thermostat_type
      LOGICAL :: apply_general_thermo, apply_thermo_adiabatic, apply_thermo_baro, &
         apply_thermo_shell, explicit_adiabatic_fast, explicit_adiabatic_slow, explicit_baro, &
         explicit_barostat_section, explicit_part, explicit_shell, save_mem, shell_adiabatic, &
         shell_present
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qmmm_env_type), POINTER                       :: qmmm_env
      TYPE(section_vals_type), POINTER :: adiabatic_fast_section, adiabatic_slow_section, &
         barostat_section, print_section, region_section_fast, region_section_slow, &
         region_sections, thermo_baro_section, thermo_part_section, thermo_shell_section, &
         work_section

      NULLIFY (qmmm_env, cell)
      ALLOCATE (thermostats)
      CALL allocate_thermostats(thermostats)
      adiabatic_fast_section => section_vals_get_subs_vals(md_section, "ADIABATIC_DYNAMICS%THERMOSTAT_FAST")
      adiabatic_slow_section => section_vals_get_subs_vals(md_section, "ADIABATIC_DYNAMICS%THERMOSTAT_SLOW")
      thermo_part_section => section_vals_get_subs_vals(md_section, "THERMOSTAT")
      thermo_shell_section => section_vals_get_subs_vals(md_section, "SHELL%THERMOSTAT")
      thermo_baro_section => section_vals_get_subs_vals(md_section, "BAROSTAT%THERMOSTAT")
      barostat_section => section_vals_get_subs_vals(md_section, "BAROSTAT")
      print_section => section_vals_get_subs_vals(md_section, "PRINT")

      CALL force_env_get(force_env, qmmm_env=qmmm_env, subsys=subsys, cell=cell)
      CALL section_vals_get(barostat_section, explicit=explicit_barostat_section)
      CALL section_vals_val_get(global_section, "SAVE_MEM", l_val=save_mem)
      CALL section_vals_get(thermo_part_section, explicit=explicit_part)
      CALL section_vals_get(thermo_shell_section, explicit=explicit_shell)
      CALL section_vals_get(thermo_baro_section, explicit=explicit_baro)
      CALL section_vals_get(adiabatic_fast_section, explicit=explicit_adiabatic_fast)
      CALL section_vals_get(adiabatic_slow_section, explicit=explicit_adiabatic_slow)

      apply_thermo_adiabatic = (simpar%ensemble == nvt_adiabatic_ensemble)

      apply_thermo_baro = (simpar%ensemble == npt_f_ensemble) .OR. &
                          (simpar%ensemble == npt_ia_ensemble) .OR. &
                          (simpar%ensemble == npt_i_ensemble) .AND. &
                          (.NOT. apply_thermo_adiabatic)

      apply_general_thermo = apply_thermo_baro .OR. (simpar%ensemble == nvt_ensemble) .AND. &
                             (.NOT. apply_thermo_adiabatic)

      apply_thermo_shell = (simpar%ensemble == nve_ensemble) .OR. &
                           (simpar%ensemble == nvt_ensemble) .OR. &
                           (simpar%ensemble == npt_f_ensemble) .OR. &
                           (simpar%ensemble == npt_i_ensemble) .OR. &
                           (simpar%ensemble == npt_ia_ensemble) .OR. &
                           (simpar%ensemble == npe_i_ensemble) .OR. &
                           (simpar%ensemble == npe_f_ensemble) .AND. &
                           (.NOT. apply_thermo_adiabatic)

      binary_restart_file_name = ""
      CALL section_vals_val_get(force_env%root_section, "EXT_RESTART%BINARY_RESTART_FILE_NAME", &
                                c_val=binary_restart_file_name)

      ! Compute Degrees of Freedom
      IF (simpar%ensemble == nvt_adiabatic_ensemble) THEN
         CALL cite_reference(VandeVondele2002)
         region = do_region_global
         region_section_fast => section_vals_get_subs_vals(adiabatic_fast_section, "DEFINE_REGION")
         region_section_slow => section_vals_get_subs_vals(adiabatic_slow_section, "DEFINE_REGION")
         IF (explicit_adiabatic_fast) CALL section_vals_val_get(adiabatic_fast_section, "REGION", i_val=region)
         IF (explicit_adiabatic_slow) CALL section_vals_val_get(adiabatic_slow_section, "REGION", i_val=region)
         CALL cp_subsys_get(subsys, molecule_kinds=molecule_kinds, local_molecules=local_molecules, &
                            molecules=molecules, gci=gci, particles=particles)
         CALL compute_nfree(cell, simpar, molecule_kinds%els, &
                            print_section, particles, gci)
         IF (explicit_adiabatic_fast .AND. explicit_adiabatic_slow) THEN
            IF (apply_thermo_adiabatic) THEN
               ALLOCATE (thermostats%thermostat_fast)
               CALL create_thermostat_type(thermostats%thermostat_fast, simpar, adiabatic_fast_section, &
                                           label="FAST")
               ALLOCATE (thermostats%thermostat_slow)
               CALL create_thermostat_type(thermostats%thermostat_slow, simpar, adiabatic_slow_section, &
                                           label="SLOW")
               CALL setup_adiabatic_thermostat_info(thermostats%thermostat_info_fast, &
                                                    molecule_kinds%els, local_molecules, molecules, particles, &
                                                    region, simpar%ensemble, region_sections=region_section_fast, &
                                                    qmmm_env=qmmm_env)

               CALL setup_adiabatic_thermostat_info(thermostats%thermostat_info_slow, &
                                                    molecule_kinds%els, local_molecules, molecules, particles, &
                                                    region, simpar%ensemble, region_sections=region_section_slow, &
                                                    qmmm_env=qmmm_env)

               ! Initialize or possibly restart Nose on Particles
               work_section => section_vals_get_subs_vals(adiabatic_fast_section, "NOSE")
               CALL initialize_nhc_fast(thermostats%thermostat_info_fast, simpar, local_molecules, &
                                        molecules%els, molecule_kinds%els, para_env, globenv, &
                                        thermostats%thermostat_fast%nhc, nose_section=work_section, gci=gci, &
                                        save_mem=save_mem)
               work_section => section_vals_get_subs_vals(adiabatic_slow_section, "NOSE")
               CALL initialize_nhc_slow(thermostats%thermostat_info_slow, simpar, local_molecules, &
                                        molecules%els, molecule_kinds%els, para_env, globenv, &
                                        thermostats%thermostat_slow%nhc, nose_section=work_section, gci=gci, &
                                        save_mem=save_mem)
            END IF
         ELSE
            CALL cp_warn(__LOCATION__, &
                         "Adiabatic Thermostat has been defined but the ensemble provided "// &
                         "does not support thermostat for Particles! Ignoring thermostat input.")
         END IF
         CALL release_thermostat_info(thermostats%thermostat_info_fast)
         DEALLOCATE (thermostats%thermostat_info_fast)
         CALL release_thermostat_info(thermostats%thermostat_info_slow)
         DEALLOCATE (thermostats%thermostat_info_fast)
      ELSE
         region = do_region_global
         region_sections => section_vals_get_subs_vals(thermo_part_section, "DEFINE_REGION")
         IF (explicit_part) CALL section_vals_val_get(thermo_part_section, "REGION", i_val=region)
         CALL cp_subsys_get(subsys, molecule_kinds=molecule_kinds, local_molecules=local_molecules, &
                            molecules=molecules, gci=gci, particles=particles)
         CALL compute_degrees_of_freedom(thermostats, cell, simpar, molecule_kinds%els, &
                                         local_molecules, molecules, particles, print_section, region_sections, gci, &
                                         region, qmmm_env)

         ! Particles
         ! For constant temperature ensembles the thermostat is activated by default
         IF (explicit_part) THEN
            IF (apply_general_thermo) THEN
               ALLOCATE (thermostats%thermostat_part)
               CALL create_thermostat_type(thermostats%thermostat_part, simpar, thermo_part_section, &
                                           label="PARTICLES")
               ! Initialize thermostat
               IF (thermostats%thermostat_part%type_of_thermostat == do_thermo_nose) THEN
                  ! Initialize or possibly restart Nose on Particles
                  work_section => section_vals_get_subs_vals(thermo_part_section, "NOSE")
                  CALL initialize_nhc_part(thermostats%thermostat_info_part, simpar, local_molecules, &
                                           molecules%els, molecule_kinds%els, para_env, globenv, &
                                           thermostats%thermostat_part%nhc, nose_section=work_section, gci=gci, &
                                           save_mem=save_mem, binary_restart_file_name=binary_restart_file_name)
               ELSE IF (thermostats%thermostat_part%type_of_thermostat == do_thermo_csvr) THEN
                  ! Initialize or possibly restart CSVR thermostat on Particles
                  work_section => section_vals_get_subs_vals(thermo_part_section, "CSVR")
                  CALL initialize_csvr_part(thermostats%thermostat_info_part, simpar, local_molecules, &
                                            molecules%els, molecule_kinds%els, para_env, &
                                            thermostats%thermostat_part%csvr, csvr_section=work_section, &
                                            gci=gci)
               ELSE IF (thermostats%thermostat_part%type_of_thermostat == do_thermo_al) THEN
                  ! Initialize or possibly restart Ad-Langevin thermostat on Particles
                  work_section => section_vals_get_subs_vals(thermo_part_section, "AD_LANGEVIN")
                  CALL initialize_al_part(thermostats%thermostat_info_part, simpar, local_molecules, &
                                          molecules%els, molecule_kinds%els, para_env, &
                                          thermostats%thermostat_part%al, al_section=work_section, &
                                          gci=gci)
               ELSE IF (thermostats%thermostat_part%type_of_thermostat == do_thermo_gle) THEN
                  ! Initialize or possibly restart GLE thermostat on Particles
                  work_section => section_vals_get_subs_vals(thermo_part_section, "GLE")
                  CALL initialize_gle_part(thermostats%thermostat_info_part, simpar, local_molecules, &
                                           molecules%els, molecule_kinds%els, para_env, &
                                           thermostats%thermostat_part%gle, gle_section=work_section, &
                                           gci=gci, save_mem=save_mem)
               END IF
               CALL thermostat_info(thermostats%thermostat_part, "PARTICLES", thermo_part_section, &
                                    simpar, para_env)
            ELSE
               CALL cp_warn(__LOCATION__, &
                            "Thermostat for Particles has been defined but the ensemble provided "// &
                            "does not support thermostat for Particles! Ignoring thermostat input.")
            END IF
         ELSE IF (apply_general_thermo) THEN
            CALL cp_abort(__LOCATION__, &
                          "One constant temperature ensemble has been required, but no thermostat for the "// &
                          "particles has been defined. You may want to change your input and add a "// &
                          "THERMOSTAT section in the MD section.")
         END IF

         ! Core-Shell Model
         CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els, shell_present=shell_present)
         IF (shell_present) THEN
            IF (explicit_shell) THEN
               ! The thermostat is activated only if explicitely required
               ! It can be used to thermalize the shell-core motion when the temperature is not constant (nve, npe)
               IF (apply_thermo_shell) THEN
                  ALLOCATE (thermostats%thermostat_shell)
                  CALL create_thermostat_type(thermostats%thermostat_shell, simpar, thermo_shell_section, &
                                              label="SHELL")
                  CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els, shell_adiabatic=shell_adiabatic)
                  region_sections => section_vals_get_subs_vals(thermo_shell_section, "DEFINE_REGION")
                  CALL section_vals_val_get(thermo_shell_section, "REGION", i_val=region)
                  CALL setup_thermostat_info( &
                     thermostats%thermostat_info_shell, molecule_kinds%els, &
                     local_molecules, molecules, particles, region, simpar%ensemble, shell=shell_adiabatic, &
                     region_sections=region_sections, qmmm_env=qmmm_env)
                  IF (shell_adiabatic) THEN
                     ! Initialize thermostat
                     IF (thermostats%thermostat_shell%type_of_thermostat == do_thermo_nose) THEN
                        ! Initialize or possibly restart Nose on Shells
                        work_section => section_vals_get_subs_vals(thermo_shell_section, "NOSE")
                        CALL initialize_nhc_shell(thermostats%thermostat_info_shell, simpar, local_molecules, &
                                                  molecules%els, molecule_kinds%els, para_env, globenv, &
                                                  thermostats%thermostat_shell%nhc, nose_section=work_section, gci=gci, &
                                                  save_mem=save_mem, binary_restart_file_name=binary_restart_file_name)
                     ELSE IF (thermostats%thermostat_shell%type_of_thermostat == do_thermo_csvr) THEN
                        ! Initialize or possibly restart CSVR thermostat on Shells
                        work_section => section_vals_get_subs_vals(thermo_shell_section, "CSVR")
                        CALL initialize_csvr_shell(thermostats%thermostat_info_shell, simpar, local_molecules, &
                                                   molecules%els, molecule_kinds%els, para_env, &
                                                   thermostats%thermostat_shell%csvr, csvr_section=work_section, gci=gci)
                     END IF
                     CALL thermostat_info(thermostats%thermostat_shell, "CORE-SHELL", thermo_shell_section, &
                                          simpar, para_env)
                  ELSE
                     CALL cp_warn(__LOCATION__, &
                                  "Thermostat for Core-Shell motion only with adiabatic shell-model. "// &
                                  "Continuing calculation ignoring the thermostat info! No Thermostat "// &
                                  "applied to Shells!")
                     CALL release_thermostat_type(thermostats%thermostat_shell)
                     DEALLOCATE (thermostats%thermostat_shell)
                     CALL release_thermostat_info(thermostats%thermostat_info_shell)
                     DEALLOCATE (thermostats%thermostat_info_shell)
                  END IF
               ELSE
                  CALL cp_warn(__LOCATION__, &
                               "Thermostat for Shells has been defined but for the selected ensemble the adiabatic "// &
                               "shell model has not been implemented! Ignoring thermostat input.")
               END IF
            END IF
         ELSE IF (explicit_shell) THEN
            CALL cp_warn(__LOCATION__, &
                         "Thermostat for Shells has been defined but the system provided "// &
                         "does not contain any Shells! Ignoring thermostat input.")
         END IF

         ! Barostat Temperature (not necessarily to be controlled by a thermostat)
         IF (explicit_barostat_section) THEN
            simpar%temp_baro_ext = simpar%temp_ext
            CALL section_vals_val_get(md_section, "BAROSTAT%TEMPERATURE", n_rep_val=n_rep)
            IF (n_rep /= 0) THEN
               CALL section_vals_val_get(md_section, "BAROSTAT%TEMPERATURE", r_val=simpar%temp_baro_ext)
               CPASSERT(simpar%temp_baro_ext >= 0.0_dp)
            END IF

            ! Setup Barostat Thermostat
            IF (apply_thermo_baro) THEN
               ! Check if we use the same thermostat as particles
               CALL section_vals_val_get(thermo_baro_section, "TYPE", i_val=thermostat_type)
               work_section => thermo_baro_section
               IF (thermostat_type == do_thermo_same_as_part) work_section => thermo_part_section

               ALLOCATE (thermostats%thermostat_baro)
               CALL create_thermostat_type(thermostats%thermostat_baro, simpar, work_section, skip_region=.TRUE., &
                                           label="BAROSTAT")
               ! Initialize thermostat
               IF (thermostats%thermostat_baro%type_of_thermostat == do_thermo_nose) THEN
                  ! Initialize or possibly restart Nose on Barostat
                  work_section => section_vals_get_subs_vals(thermo_baro_section, "NOSE")
                  CALL initialize_nhc_baro(simpar, para_env, globenv, thermostats%thermostat_baro%nhc, &
                                           nose_section=work_section, save_mem=save_mem)
               ELSE IF (thermostats%thermostat_baro%type_of_thermostat == do_thermo_csvr) THEN
                  ! Initialize or possibly restart CSVR thermostat on Barostat
                  work_section => section_vals_get_subs_vals(thermo_baro_section, "CSVR")
                  CALL initialize_csvr_baro(simpar, thermostats%thermostat_baro%csvr, &
                                            csvr_section=work_section)
               END IF
               CALL thermostat_info(thermostats%thermostat_baro, "BAROSTAT", thermo_baro_section, &
                                    simpar, para_env)
               ! If thermostat for barostat uses a diffent kind than the one of the particles
               ! let's update infos in the input structure..
               IF (thermostat_type == do_thermo_same_as_part) THEN
                  CALL update_thermo_baro_section(thermostats%thermostat_baro, thermo_baro_section)
               END IF
            ELSE
               IF (explicit_baro) THEN
                  CALL cp_warn(__LOCATION__, &
                               "Thermostat for Barostat has been defined but the ensemble provided "// &
                               "does not support thermostat for Barostat! Ignoring thermostat input.")
               END IF
               ! Let's remove the section
               CALL section_vals_remove_values(thermo_baro_section)
            END IF
         END IF

         ! Release the thermostats info..
         CALL release_thermostat_info(thermostats%thermostat_info_part)
         DEALLOCATE (thermostats%thermostat_info_part)
         CALL release_thermostat_info(thermostats%thermostat_info_shell)
         DEALLOCATE (thermostats%thermostat_info_shell)

      END IF ! Adiabitic_NVT screening
      ! If no thermostats have been allocated deallocate the full structure
      IF ((.NOT. ASSOCIATED(thermostats%thermostat_part)) .AND. &
          (.NOT. ASSOCIATED(thermostats%thermostat_shell)) .AND. &
          (.NOT. ASSOCIATED(thermostats%thermostat_baro)) .AND. &
          (.NOT. ASSOCIATED(thermostats%thermostat_fast)) .AND. &
          (.NOT. ASSOCIATED(thermostats%thermostat_slow))) THEN
         CALL release_thermostats(thermostats)
         DEALLOCATE (thermostats)
      END IF

   END SUBROUTINE create_thermostats

! **************************************************************************************************
!> \brief ...
!> \param thermostat ...
!> \param section ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE update_thermo_baro_section(thermostat, section)
      TYPE(thermostat_type), POINTER                     :: thermostat
      TYPE(section_vals_type), POINTER                   :: section

      TYPE(section_vals_type), POINTER                   :: work_section

      CALL section_vals_val_set(section, "TYPE", i_val=thermostat%type_of_thermostat)
      SELECT CASE (thermostat%type_of_thermostat)
      CASE (do_thermo_nose)
         work_section => section_vals_get_subs_vals(section, "NOSE")
         CALL section_vals_val_set(work_section, "LENGTH", i_val=thermostat%nhc%nhc_len)
         CALL section_vals_val_set(work_section, "YOSHIDA", i_val=thermostat%nhc%nyosh)
         CALL section_vals_val_set(work_section, "TIMECON", r_val=thermostat%nhc%tau_nhc)
         CALL section_vals_val_set(work_section, "MTS", i_val=thermostat%nhc%nc)
      CASE (do_thermo_csvr)
         work_section => section_vals_get_subs_vals(section, "CSVR")
         CALL section_vals_val_set(work_section, "TIMECON", r_val=thermostat%csvr%tau_csvr)
      CASE (do_thermo_al)
         work_section => section_vals_get_subs_vals(section, "AD_LANGEVIN")
         CALL section_vals_val_set(work_section, "TIMECON_NH", r_val=thermostat%al%tau_nh)
         CALL section_vals_val_set(work_section, "TIMECON_LANGEVIN", r_val=thermostat%al%tau_langevin)
      END SELECT
   END SUBROUTINE update_thermo_baro_section

! **************************************************************************************************
!> \brief ...
!> \param thermostat ...
!> \param label ...
!> \param section ...
!> \param simpar ...
!> \param para_env ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE thermostat_info(thermostat, label, section, simpar, para_env)
      TYPE(thermostat_type), POINTER                     :: thermostat
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: iw
      REAL(KIND=dp)                                      :: kin_energy, pot_energy, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, section, "PRINT%THERMOSTAT_INFO", extension=".log")
      ! Total Tehrmostat Energy
      CALL get_thermostat_energies(thermostat, pot_energy, kin_energy, para_env)
      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            'THERMOSTAT| Thermostat information for '//TRIM(label)
         SELECT CASE (thermostat%type_of_thermostat)
         CASE (do_thermo_nose)
            WRITE (iw, '(T2,A,T63,A)') &
               'THERMOSTAT| Type of thermostat', 'Nose-Hoover-Chains'
            WRITE (iw, '(T2,A,T71,I10)') &
               'THERMOSTAT| Nose-Hoover-Chain length', thermostat%nhc%nhc_len
            tmp = cp_unit_from_cp2k(thermostat%nhc%tau_nhc, 'fs')
            WRITE (iw, '(T2,A,T61,F20.6)') &
               'THERMOSTAT| Nose-Hoover-Chain time constant [fs]', tmp
            WRITE (iw, '(T2,A,T71,I10)') &
               'THERMOSTAT| Order of Yoshida integrator', thermostat%nhc%nyosh
            WRITE (iw, '(T2,A,T71,I10)') &
               'THERMOSTAT| Number of multiple time steps', thermostat%nhc%nc
            WRITE (iw, '(T2,A,T61,E20.12)') &
               'THERMOSTAT| Initial potential energy', pot_energy, &
               'THERMOSTAT| Initial kinetic energy', kin_energy
         CASE (do_thermo_csvr)
            WRITE (iw, '(T2,A,T44,A)') &
               'THERMOSTAT| Type of thermostat', 'Canonical Sampling/Velocity Rescaling'
            tmp = cp_unit_from_cp2k(thermostat%csvr%tau_csvr, 'fs')*0.5_dp*simpar%dt
            WRITE (iw, '(T2,A,T61,F20.6)') &
               'THERMOSTAT| CSVR time constant [fs]', tmp
            WRITE (iw, '(T2,A,T61,E20.12)') &
               'THERMOSTAT| Initial kinetic energy', kin_energy
         CASE (do_thermo_al)
            WRITE (iw, '(T2,A,T44,A)') &
               'THERMOSTAT| Type of thermostat', 'Adaptive Langevin'
            tmp = cp_unit_from_cp2k(thermostat%al%tau_nh, 'fs')
            WRITE (iw, '(T2,A,T61,F20.6)') &
               'THERMOSTAT| Time constant of Nose-Hoover part [fs]', tmp
            tmp = cp_unit_from_cp2k(thermostat%al%tau_langevin, 'fs')
            WRITE (iw, '(T2,A,T61,F20.6)') &
               'THERMOSTAT| Time constant of Langevin part [fs]', tmp
         END SELECT
         WRITE (iw, '(T2,A)') &
            'THERMOSTAT| End of thermostat information for '//TRIM(label)
      END IF
      CALL cp_print_key_finished_output(iw, logger, section, "PRINT%THERMOSTAT_INFO")

   END SUBROUTINE thermostat_info

! **************************************************************************************************
!> \brief ...
!> \param thermostat ...
!> \param npt ...
!> \param group ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE apply_thermostat_baro(thermostat, npt, group)
      TYPE(thermostat_type), POINTER                     :: thermostat
      TYPE(npt_info_type), DIMENSION(:, :), POINTER      :: npt

      CLASS(mp_comm_type), INTENT(IN)                     :: group

      IF (ASSOCIATED(thermostat)) THEN
         IF (thermostat%type_of_thermostat == do_thermo_nose) THEN
            ! Apply Nose-Hoover Thermostat
            CPASSERT(ASSOCIATED(thermostat%nhc))
            CALL lnhc_barostat(thermostat%nhc, npt, group)
         ELSE IF (thermostat%type_of_thermostat == do_thermo_csvr) THEN
            ! Apply CSVR Thermostat
            CPASSERT(ASSOCIATED(thermostat%csvr))
            CALL csvr_barostat(thermostat%csvr, npt, group)
         END IF
      END IF
   END SUBROUTINE apply_thermostat_baro

! **************************************************************************************************
!> \brief ...
!> \param thermostat ...
!> \param force_env ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param local_particles ...
!> \param group ...
!> \param shell_adiabatic ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE apply_thermostat_particles(thermostat, force_env, molecule_kind_set, molecule_set, &
                                         particle_set, local_molecules, local_particles, &
                                         group, shell_adiabatic, shell_particle_set, &
                                         core_particle_set, vel, shell_vel, core_vel)

      TYPE(thermostat_type), POINTER                     :: thermostat
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles

      CLASS(mp_comm_type), INTENT(IN)                     :: group
      LOGICAL, INTENT(IN), OPTIONAL                      :: shell_adiabatic
      TYPE(particle_type), OPTIONAL, POINTER             :: shell_particle_set(:), &
                                                            core_particle_set(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: vel(:, :), shell_vel(:, :), &
                                                            core_vel(:, :)

      IF (ASSOCIATED(thermostat)) THEN
         IF (thermostat%type_of_thermostat == do_thermo_nose) THEN
            ! Apply Nose-Hoover Thermostat
            CPASSERT(ASSOCIATED(thermostat%nhc))
            CALL lnhc_particles(thermostat%nhc, molecule_kind_set, molecule_set, &
                                particle_set, local_molecules, group, shell_adiabatic, shell_particle_set, &
                                core_particle_set, vel, shell_vel, core_vel)
         ELSE IF (thermostat%type_of_thermostat == do_thermo_csvr) THEN
            ! Apply CSVR Thermostat
            CPASSERT(ASSOCIATED(thermostat%csvr))
            CALL csvr_particles(thermostat%csvr, molecule_kind_set, molecule_set, &
                                particle_set, local_molecules, group, shell_adiabatic, shell_particle_set, &
                                core_particle_set, vel, shell_vel, core_vel)
         ELSE IF (thermostat%type_of_thermostat == do_thermo_al) THEN
            ! Apply AD_LANGEVIN Thermostat
            CPASSERT(ASSOCIATED(thermostat%al))
            CALL al_particles(thermostat%al, force_env, molecule_kind_set, molecule_set, &
                              particle_set, local_molecules, local_particles, group, vel)
         ELSE IF (thermostat%type_of_thermostat == do_thermo_gle) THEN
            ! Apply GLE Thermostat
            CPASSERT(ASSOCIATED(thermostat%gle))
            CALL gle_particles(thermostat%gle, molecule_kind_set, molecule_set, &
                               particle_set, local_molecules, group, shell_adiabatic, shell_particle_set, &
                               core_particle_set, vel, shell_vel, core_vel)
         END IF
      END IF
   END SUBROUTINE apply_thermostat_particles

! **************************************************************************************************
!> \brief ...
!> \param thermostat ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param group ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE apply_thermostat_shells(thermostat, atomic_kind_set, particle_set, &
                                      local_particles, group, shell_particle_set, core_particle_set, vel, shell_vel, &
                                      core_vel)

      TYPE(thermostat_type), POINTER                     :: thermostat
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles

      CLASS(mp_comm_type), INTENT(IN)                     :: group
      TYPE(particle_type), OPTIONAL, POINTER             :: shell_particle_set(:), &
                                                            core_particle_set(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: vel(:, :), shell_vel(:, :), &
                                                            core_vel(:, :)

      IF (ASSOCIATED(thermostat)) THEN
         IF (thermostat%type_of_thermostat == do_thermo_nose) THEN
            ! Apply Nose-Hoover Thermostat
            CPASSERT(ASSOCIATED(thermostat%nhc))
            CALL lnhc_shells(thermostat%nhc, atomic_kind_set, particle_set, local_particles, &
                             group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel)
         ELSE IF (thermostat%type_of_thermostat == do_thermo_csvr) THEN
            ! Apply CSVR Thermostat
            CPASSERT(ASSOCIATED(thermostat%csvr))
            CALL csvr_shells(thermostat%csvr, atomic_kind_set, particle_set, local_particles, &
                             group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel)
         END IF
      END IF
   END SUBROUTINE apply_thermostat_shells

END MODULE thermostat_methods
