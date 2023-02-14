!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Nonequilibrium structure: module containing nonequilibrium forces available for MD
!> \author James Reid [JR], Emily Kahl [EVK] - University of Queensland - 2023
! *****************************************************************************
MODULE nonequilibrium_types
  USE input_constants,                 ONLY: do_ne_colour,&
                                              do_ne_field
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE simpar_types,                    ONLY: simpar_type
  USE particle_list_types,             ONLY: particle_list_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
#include "../base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: NE_type,&
!            allocate_NE,&
!            release_NE,&
!            retain_NE,&
            create_NE_type
!            release_NE_type,&
!            NE_info_type,&
!            release_NE_info,&
!            set_NE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'Non-equilibrium_methods'

  ! Single NE_type JR - should trim down unused parts!!
! *****************************************************************************
  TYPE NE_type
     INTEGER :: id_nr, ref_count, type_of_NE
     CHARACTER(LEN=default_string_length) :: label
     REAL(KIND=dp)                            :: FX, FY, FZ
      TYPE(section_vals_type), POINTER :: section
  END TYPE NE_type

CONTAINS

   SUBROUTINE create_NE_type(nevar, simpar, md_section, force_env)
      TYPE(NE_type), POINTER                             :: nevar
      TYPE(section_vals_type), POINTER                   :: md_section, section
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(cp_subsys_type), POINTER                      :: subsys

      LOGICAL                                            :: explicit, shell_present
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_list_type), POINTER                  :: particles, shell_particles
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      INTEGER                                            :: i, ipart

      NULLIFY (subsys)
      NULLIFY (atomic_kinds, atomic_kind_set, particles, shell_particles)
      simpar%do_NEforce = .FALSE.
      ! TODO EVK: This check is completely busted. The CPASSERT is getting flagged
      ! by Valgrind as depending on uninitialsed memory, despite "check"
      ! being assigned to right here. But if I throw in the print
      ! statement then everything's hunky dory? This is wack.
      !check = .NOT. ASSOCIATED(nevar)
      !print*, "================================= Check = ", check
      !CPASSERT(check)
      IF (ASSOCIATED(nevar)) THEN
        CPABORT("NE_type already allocated")
      END IF

      section => section_vals_get_subs_vals(md_section, "NE")
      CALL section_vals_get(section, explicit=explicit)
      IF (explicit) THEN 
        ALLOCATE(nevar)
        nevar%id_nr = 1
        nevar%ref_count = 1
        simpar%do_NEforce = .TRUE.
        nevar%section => section

        CALL section_vals_val_get(section, "NETYPE", i_val=nevar%type_of_NE)
        CALL section_vals_val_get(section, "FX", r_val=nevar%FX)
        CALL section_vals_val_get(section, "FY", r_val=nevar%FY)
        CALL section_vals_val_get(section, "FZ", r_val=nevar%FZ)

        ! Now get the list of all particles, via the force subsys
        CALL force_env_get(force_env, subsys=subsys)
        CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles, &
                           shell_particles=shell_particles)
        atomic_kind_set => atomic_kinds%els

        ! Colour charge is not tested for core-shell model, so bail out
        ! if shell_present
        CALL get_atomic_kind_set(atomic_kind_set, shell_present=shell_present)
        IF (shell_present) THEN
          ! TODO EVK: should this be CPABORT or CPASSERT?
          CPABORT("NETYPE not compatible with shells: particles only")
        END IF

        ! Assign colour charges to particles
        particles%els(:)%ccharge = 0
        CALL section_vals_val_get(section, "NELIST", i_vals=tmplist)
        IF (nevar%type_of_NE .EQ. do_ne_colour .AND. mod(SIZE(tmplist), 2) .NE. 0) THEN
          CPABORT("NELIST: colour field must have an even number of particles")
        END IF

        DO i = 1, SIZE(tmplist)
          ipart = tmplist(i)
          IF (nevar%type_of_NE .EQ. do_ne_colour) THEN
            particles%els(ipart)%ccharge = (-1)**(i)
          ELSE IF (nevar%type_of_NE .EQ. do_ne_field) THEN
            particles%els(ipart)%ccharge = 1
          END IF
        END DO

        ! Check the colour charges all sum to zero. This should never happen
        if (nevar%type_of_NE .EQ. do_ne_colour .AND. sum(particles%els(:)%ccharge) .NE. 0) THEN
          CPABORT("NELIST: colour charges do not sum to zero")
        END IF
        
        CALL NE_info(nevar, particles)
      END IF

   END SUBROUTINE create_NE_type

   SUBROUTINE NE_info(nevar, particles)
      TYPE(NE_type), POINTER                             :: nevar
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: i, iw
      LOGICAL                                            :: new_file

      NULLIFY (logger)
      logger => cp_get_default_logger()

      IF (ASSOCIATED(nevar)) THEN
        iw = cp_print_key_unit_nr(logger, nevar%section, "PRINT%NE_INFO", &
                                 extension=".log", is_new_file=new_file)

        IF (iw > 0) THEN
          IF (new_file) THEN
            WRITE ( iw, '( /,A )' ) ' NON-EQUILIBRIUM| NE Info for MD'
          END IF

          IF (nevar%type_of_NE.eq.do_ne_colour) THEN
            WRITE ( iw, '( A,T69,A )' ) ' NON-EQUILIBRIUM| Type of thermostat','Colour Field'
          ELSE
            WRITE ( iw, '( A,T65,A )' ) ' NON-EQUILIBRIUM| Type of thermostat','Poiseuille Field'
          END IF

          WRITE ( iw, '( A,A5,A,T69,F12.9 )' ) &
               ' NON-EQUILIBRIUM| Force in X Direction [', 'au_f' ,'] ',nevar%FX
          WRITE ( iw, '( A,A5,A,T69,F12.9 )' ) &
               ' NON-EQUILIBRIUM| Force in Y Direction [', 'au_f' ,'] ',nevar%FY
          WRITE ( iw, '( A,A5,A,T69,F12.9 )' ) &
               ' NON-EQUILIBRIUM| Force in Z Direction [', 'au_f' ,'] ',nevar%FZ

          WRITE ( iw, '(A,/)' ) ' NON-EQUILIBRIUM| NE "Charges":'

          WRITE ( iw, '( A,T72,A,/ )' ) ' Atom', 'NE charge'
          DO i=1,size(particles%els(:))
            WRITE ( iw, '( I6,T79,I2 )' ) i, particles%els(i)%ccharge
          END DO
          WRITE(iw,'(/,/)')
        END IF
        
        CALL cp_print_key_finished_output(iw, logger, nevar%section,"PRINT%NE_INFO")

      END IF

   END SUBROUTINE NE_info

END MODULE nonequilibrium_types
