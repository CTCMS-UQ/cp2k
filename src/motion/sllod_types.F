MODULE sllod_types
  !USE input_constants,                 ONLY:!, do_XY!, different types of
  !shear force
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
  PUBLIC :: SLLOD_type, &
            create_SLLOD_type, &
            SLLOD_info

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'Non-equilibrium_methods' !could this be sllod_methods?


  ! Adapted from JR and EK nonequilibrium_types.F for SLLOD - LM
! *****************************************************************************
  TYPE SLLOD_type
     INTEGER :: id_nr, ref_count, type_of_SLLOD, remap_prd
     CHARACTER(LEN=default_string_length) :: label
     REAL(KIND=dp)                            :: XY
     TYPE(section_vals_type), POINTER :: section
  END TYPE SLLOD_type

CONTAINS

   SUBROUTINE create_SLLOD_type(sllodvar, simpar, md_section, force_env)
      TYPE(SLLOD_type), POINTER                          :: sllodvar
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
      simpar%do_SLLOD = .FALSE.
      ! TODO EVK: CPASSERT issue as in nonequilibrium_types.F
      IF (ASSOCIATED(sllodvar)) THEN
        CPABORT("SLLOD_type already allocated")
      END IF

      section => section_vals_get_subs_vals(md_section, "SLLOD")
      CALL section_vals_get(section, explicit=explicit)
      IF (explicit) THEN 
        ALLOCATE(sllodvar)
        sllodvar%id_nr = 1
        sllodvar%ref_count = 1
        simpar%do_SLLOD = .TRUE.
        sllodvar%section => section

        !CALL section_vals_val_get(section, "SHEARTYPE", i_val=sllodvar%type_of_SLLOD)
        ! or just directly write in gradu in input?
        CALL section_vals_val_get(section, "XY", r_val=sllodvar%XY)
        CALL section_vals_val_get(section, "REMAP_PERIOD", i_val=sllodvar%remap_prd)
        ! Now get the list of all particles, via the force subsys
        CALL force_env_get(force_env, subsys=subsys)
        CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles, &
                           shell_particles=shell_particles)
        atomic_kind_set => atomic_kinds%els

        !bail out if core-shell model is present since not tested
        CALL get_atomic_kind_set(atomic_kind_set, shell_present=shell_present)
       IF (shell_present) THEN
         CPABORT("SLLOD not compatible with shells: particles only")
       END IF

        ! Turn particles on or off for sllod
        CALL section_vals_val_get(section, "SLLODLIST", i_vals=tmplist)

        DO i = 1, SIZE(tmplist)
          ipart = tmplist(i)
          particles%els(ipart)%sllodswitch = .TRUE.
        END DO

        
        CALL SLLOD_info(sllodvar, particles)
       END IF
   END SUBROUTINE create_SLLOD_type

   SUBROUTINE SLLOD_info(sllodvar, particles)
      TYPE(SLLOD_type), POINTER                          :: sllodvar
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: i, iw
      LOGICAL                                            :: new_file

      NULLIFY (logger)
      logger => cp_get_default_logger()

      IF (ASSOCIATED(sllodvar)) THEN
        iw = cp_print_key_unit_nr(logger, sllodvar%section, "PRINT%SLLOD_INFO", &
                                 extension=".log", is_new_file=new_file)

        IF (iw > 0) THEN
          IF (new_file) THEN
            WRITE ( iw, '( /,A )' ) ' NON-EQUILIBRIUM| Info for SLLOD MD'
          END IF

!TODO - DIFFERENT TYPES OF SHEAR IN DIFFERENT PLANES 
!          IF (sllodvar%type_of_SHEAR.eq.XY??) THEN
!            WRITE ( iw, '( A,T69,A )' ) ' NON-EQUILIBRIUM| Type of strain','Colour Field'
!          ELSE IF 
!            WRITE ( iw, '( A,T65,A )' ) ' NON-EQUILIBRIUM| Type of thermostat','Poiseuille Field'
!          END IF

          WRITE ( iw, '( A,A5,A,T69,F12.9 )' ) &
               ' NON-EQUILIBRIUM| XY STRAIN RATE [', 'au_f' ,'] ',sllodvar%XY
          WRITE ( iw, '(T2,A,T71,I10)' ) &
               ' NON-EQUILIBRIUM| REMAPPING PERIOD', sllodvar%remap_prd
          WRITE ( iw, '(A,/)' ) 'SLLOD thermostated atoms:'


          WRITE ( iw, '( A,T72,A,/ )' ) ' Atom', 'SLLOD ON (1)/OFF (0)'
          DO i=1,size(particles%els(:))
            IF(particles%els(i)%sllodswitch) THEN
            WRITE ( iw, '( I6,T79,I2 )' ) i, 1
            ELSE
            WRITE ( iw, '( I6,T79,I2 )' ) i, 0
            END IF
          END DO
          WRITE(iw,'(/,/)')
        END IF
        
        CALL cp_print_key_finished_output(iw, logger, sllodvar%section,"PRINT%SLLOD_INFO")

      END IF

   END SUBROUTINE SLLOD_info
END MODULE sllod_types
