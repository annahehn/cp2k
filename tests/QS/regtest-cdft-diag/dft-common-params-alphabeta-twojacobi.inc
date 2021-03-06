&DFT
  @IF ( ${BECKE_ACTIVE} == TRUE )
    @include becke_qs_twocons_twojacobi.inc
  @ENDIF
  @IF ( ${BECKE_ACTIVE} == FALSE )
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
    &END QS
  @ENDIF
  BASIS_SET_FILE_NAME  BASIS_MOLOPT
  POTENTIAL_FILE_NAME  POTENTIAL
  @IF ( ${RESTART_WFN} == TRUE )
      WFN_RESTART_FILE_NAME ${WFN_FILE}
  @ENDIF
  LSD
  CHARGE ${CHARGE}
  &MGRID
    CUTOFF 100
    NGRIDS 5
  &END MGRID
  &POISSON
    PERIODIC NONE
    POISSON_SOLVER WAVELET
  &END POISSON
  &SCF
    @IF ( ${RESTART_WFN} == TRUE )
      SCF_GUESS RESTART
    @ENDIF
    @IF ( ${RESTART_WFN} == FALSE )
      SCF_GUESS ATOMIC
    @ENDIF
    EPS_SCF 1.0E-5
    CHOLESKY INVERSE_DBCSR
    MAX_SCF ${MAX_SCF}
    ADDED_MOS 2 2
    &DIAGONALIZATION
      ALGORITHM STANDARD
    &END DIAGONALIZATION
    &SMEAR
      METHOD FERMI_DIRAC
      ! Use extremely high electronic temperature to get fractionally occupied MOs
      ELECTRONIC_TEMPERATURE 10000
      ! CDFT couplings are only defined when the CDFT states are in the same spin state
      FIXED_MAGNETIC_MOMENT 1
    &END SMEAR
    &OUTER_SCF ON
      ! Outer_SCF section must be turned on to facilitate switching between CDFT SCF and Diag SCF loops
      EPS_SCF 1.0E-5
      MAX_SCF 0
    &END
    &PRINT
      &RESTART
        FILENAME ./${PROJECT_NAME}
        BACKUP_COPIES 0
        COMMON_ITERATION_LEVELS 1
        &EACH
          JUST_ENERGY ${WRITE_WFN}
          QS_SCF 0
        &END EACH
      &END RESTART
      &RESTART_HISTORY OFF
      &END RESTART_HISTORY
    &END PRINT
  &END SCF
  &XC
    &XC_FUNCTIONAL PBE
    &END XC_FUNCTIONAL
    &XC_GRID
      XC_DERIV SPLINE2
      XC_SMOOTH_RHO NONE
    &END XC_GRID
  &END XC
  &PRINT
    &MO
      EIGENVALUES
      OCCUPATION_NUMBERS
      COMMON_ITERATION_LEVELS 2
      ADD_LAST NUMERIC
      FILENAME ./${NAME}
      &EACH
        QS_SCF 0
        JUST_ENERGY 1
      &END EACH
    &END MO
    &MULLIKEN OFF
    &END MULLIKEN
    &HIRSHFELD OFF
    &END HIRSHFELD
    @IF ( ${WRITE_CUBE} == TRUE )
      &E_DENSITY_CUBE ON
         STRIDE 1 1 1
      &END E_DENSITY_CUBE
    @ENDIF
  &END PRINT
&END DFT
&PRINT
  &FORCES ON
  &END
&END
