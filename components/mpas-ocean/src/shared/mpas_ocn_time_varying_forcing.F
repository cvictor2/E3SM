!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing
!
!> \brief Module to handle time-varying forcing for the ocean core
!> \author Steven Brus, Xylar Asay-Davis, LANL
!> \date 2/3/2019 
!> \details
!>  Manages the wind and pressure forcing for the ocean.
!>  Based very closely on the seaice forcing module and
!>  framework written by Adrian Turner.
!
!-----------------------------------------------------------------------

module ocn_time_varying_forcing

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_forcing
  use mpas_stream_manager
  use ocn_framework_forcing
  use ocn_constants
  use ocn_config
  use ocn_diagnostics_variables
  use mpas_log, only: mpas_log_write

  implicit none

  private
  public :: &
       ocn_time_varying_forcing_init, &
       ocn_time_varying_forcing_get, &
       ocn_time_varying_forcing_write_restart_times, &
       ocn_time_varying_forcing_aot_obs_update_vel, &
       ocn_time_varying_forcing_aot_obs_update_tracers

contains

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_init
!
!> \brief Driver for ocean core forcing initialization
!> \author Steven Brus, Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_init(domain)!{{{

    type (domain_type) :: domain

    ! init the atmospheric forcing
    if (config_use_time_varying_atmospheric_forcing) then
       call init_atmospheric_forcing(domain)
    endif

    ! init other forcings, to be added later as needed
    if (config_use_time_varying_land_ice_forcing) then
      call init_land_ice_forcing(domain)
    endif

    if (config_use_time_varying_aot_forcing) then
      call init_aot_forcing(domain)
    endif

    ! init other forcings, to be added later as needed
    !if (config_use_time_varying_****_forcing) then
    !  call init_****_forcing(domain)
    !endif

  end subroutine ocn_time_varying_forcing_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing
!
!> \brief
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing(domain)!{{{

    type (domain_type) :: domain

    select case (trim(config_time_varying_atmospheric_forcing_type))
      case ("WINDPRES")
         call init_atmospheric_forcing_WINDPRES(domain)
      case default
         call mpas_log_write("Atmospheric forcing type unknown: "//trim(config_time_varying_atmospheric_forcing_type), MPAS_LOG_CRIT)
    end select

  end subroutine init_atmospheric_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing_WINDPRES
!
!> \brief Initialize the forcing objects for wind and pressure
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>  This routine calls the MPAS_forcing module subroutines that initializes
!>  the forcings type
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing_WINDPRES(domain)!{{{

    type(domain_type) :: domain

    ! create the atmospheric forcing group
    call MPAS_forcing_init_group(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_atmospheric_forcing", &            ! forcingGroupName
         domain, &                               ! domain
         config_time_varying_atmospheric_forcing_start_time, &            ! startTimeStr
         config_time_varying_atmospheric_forcing_cycle_start, &           ! forcingCycleStart
         config_time_varying_atmospheric_forcing_cycle_duration, &        ! forcingCycleDuration
         config_do_restart)!, &                  ! restart
         !.false.)                               ! forcingCycleStartInclusive

    ! windSpeedU init
    call MPAS_forcing_init_field(&
         domain % streamManager, &               ! streamManager
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_atmospheric_forcing", &            ! forcingGroupName
         "windSpeedU", &                         ! forcingName
         "atmospheric_forcing", &                ! forcingStreamID
         "timeVaryingForcing", &                 ! poolname
         "windSpeedU", &                         ! fieldname
         "linear", &                             ! interpolationType
         config_time_varying_atmospheric_forcing_reference_time, &                 ! forcingReferenceTimeStr
         config_time_varying_atmospheric_forcing_interval)!,                      ! forcingIntervalStr
         !"next")                                ! forcingInitializationType

    ! windSpeedV init
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         forcingGroupHead, &
         "ocn_atmospheric_forcing", &
         "windSpeedV", &
         "atmospheric_forcing", &
         "timeVaryingForcing", &
         "windSpeedV", &
         "linear", &
         config_time_varying_atmospheric_forcing_reference_time, &
         config_time_varying_atmospheric_forcing_interval)!, &
         !"next")

    ! atmosphericPressure init
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         forcingGroupHead, &
         "ocn_atmospheric_forcing", &
         "atmosPressure", &
         "atmospheric_forcing", &
         "timeVaryingForcing", &
         "atmosPressure", &
         "linear", &
         config_time_varying_atmospheric_forcing_reference_time, &
         config_time_varying_atmospheric_forcing_interval)!, &
         !"next")

    call MPAS_forcing_init_field_data(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_atmospheric_forcing", &            ! forcingGroupName
         domain % streamManager, &               ! streamManager
         config_do_restart, &                    ! restart
         .false.)                                ! interpolateAtInit

  end subroutine init_atmospheric_forcing_WINDPRES!}}}

!-----------------------------------------------------------------------
! runtime
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_get
!
!> \brief Retrieve forcing data during time stepping
!> \author Steven Brus, Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details
!>  This routine calls the MPAS_forcing routine that will perform the
!>  forcing data aquisition and interpolation during timestepping
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_get(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    if (config_use_time_varying_atmospheric_forcing) then

        call atmospheric_forcing(&
             streamManager, &
             domain, &
             simulationClock)

    endif

    if (config_use_time_varying_land_ice_forcing) then

        call land_ice_forcing(&
             streamManager, &
             domain, &
             simulationClock)

    endif

    if (config_use_time_varying_aot_forcing) then

        call aot_forcing(&
             streamManager, &
             domain, &
             simulationClock)

    endif

    ! call other forcings, to be added if needed
    ! if (config_use_****_forcing) then
    !   call ****_forcing(&
    !        streamManager, &
    !        domain, &
    !        simulationClock)
    ! endif


  end subroutine ocn_time_varying_forcing_get!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  atmospheric_forcing
!
!> \brief
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine atmospheric_forcing(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    type (MPAS_time_type) :: currentForcingTime

    type (block_type), pointer :: block

    type (MPAS_timeInterval_type) :: timeStepESMF

    character(len=StrKIND) :: timeStamp

    real (kind=RKIND) :: dtSim

    integer :: err

    ! convert config_dt to real
    call mpas_set_timeInterval(timeStepESMF, timeString=config_dt,ierr=err)
    call mpas_get_timeInterval(timeStepESMF, dt=dtSim)

    ! use the forcing layer to get data
    call MPAS_forcing_get_forcing(&
         forcingGroupHead, &               ! forcingGroupHead
         "ocn_atmospheric_forcing", &      ! forcingGroupName
         streamManager, &                  ! streamManager
         dtSim)                            ! dt

    call MPAS_forcing_get_forcing_time(&
         forcingGroupHead, &               ! forcingGroupHead
         "ocn_atmospheric_forcing", &      ! forcingGroupName
         currentForcingTime)               ! forcingTime
    !call mpas_get_time(curr_time=currentForcingTime, dateTimeString=timeStamp, ierr=err)
    !call mpas_log_write('Forcing time ' // trim(timeStamp))

    ! perform post forcing
    block => domain % blocklist
    do while (associated(block))

       ! apply drag law to compute wind stress from wind speed
       if (trim(config_time_varying_atmospheric_forcing_type) == "WINDPRES") then
          call post_atmospheric_forcing(block)
       endif

       block => block % next
    end do

  end subroutine atmospheric_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  post_atmospheric_forcing
!
!> \brief Calculate wind stress from wind speed and apply ramp
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>  Uses drag law from Garratt, J.R., 1977. Review of drag coefficients
!>  over oceans and continents. Monthly weather review, 105(7), pp.915-929.
!>
!
!-----------------------------------------------------------------------

  subroutine post_atmospheric_forcing(block)!{{{

    type (block_type), pointer :: block

    type (MPAS_time_type) :: currentForcingTime

    type (mpas_pool_type), pointer :: &
         mesh, &
         forcingPool, &
         timeVaryingForcingPool

    real(kind=RKIND), dimension(:), pointer :: &
         windSpeedU, &
         windSpeedV, &
         windSpeedMagnitude, &
         windStressZonal, &
         windStressMeridional, &
         atmosPressure, &
         atmosphericPressure

    real(kind=RKIND) :: &
         windStressCoefficient, &
         rhoAir, &
         ramp, &
         currentTime, &
         windStressCoefficientLimit

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)
    call MPAS_pool_get_subpool(block % structs, "timeVaryingForcing", timeVaryingForcingPool)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(timeVaryingForcingPool, "windSpeedU", windSpeedU)
    call MPAS_pool_get_array(timeVaryingForcingPool, "windSpeedV", windSpeedV)
    call MPAS_pool_get_array(timeVaryingForcingPool, "windSpeedMagnitude", windSpeedMagnitude)
    call MPAS_pool_get_array(forcingPool, "windStressZonal", windStressZonal)
    call MPAS_pool_get_array(forcingPool, "windStressMeridional", windStressMeridional)

    call MPAS_pool_get_array(timeVaryingForcingPool, "atmosPressure", atmosPressure)
    call MPAS_pool_get_array(forcingPool, "atmosphericPressure", atmosphericPressure)

    rhoAir = 1.225_RKIND
    windStressCoefficientLimit = 0.0035_RKIND

    if (daysSinceStartOfSim >= config_time_varying_atmospheric_forcing_ramp_delay) then
      ramp = tanh((2.0_RKIND*(daysSinceStartOfSim-config_time_varying_atmospheric_forcing_ramp_delay)) &
                              /config_time_varying_atmospheric_forcing_ramp)
    else
      ramp = 0.0_RKIND
    end if

    !$omp parallel
    !$omp do schedule(runtime) private(windStressCoefficient)
    do iCell = 1, nCells

       windSpeedU(iCell) = ramp*windSpeedU(iCell)
       windSpeedV(iCell) = ramp*windSpeedV(iCell)
       windSpeedMagnitude(iCell) = sqrt(windSpeedU(iCell)**2 + windSpeedV(iCell)**2)
       windStressCoefficient = (0.75_RKIND + 0.067_RKIND * windSpeedMagnitude(iCell)) * 0.001_RKIND   ! Garratt 1977
       if (windStressCoefficient > windStressCoefficientLimit) then
         windStressCoefficient = windStressCoefficientLimit
       end if

       windStressZonal(iCell) =  rhoAir * windSpeedU(iCell) * windSpeedMagnitude(iCell) * windStressCoefficient
       windStressMeridional(iCell) = rhoAir * windSpeedV(iCell) * windSpeedMagnitude(iCell) * windStressCoefficient

       atmosphericPressure(iCell) = ((1.0_RKIND-ramp)*atm_ref_pressure + ramp*atmosPressure(iCell)) - atm_ref_pressure

    enddo
    !$omp end do
    !$omp end parallel

  end subroutine post_atmospheric_forcing!}}}


!-----------------------------------------------------------------------
! restart
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_write_restart_times
!
!> \brief
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_write_restart_times(domain)!{{{

    type(domain_type) :: domain

    call MPAS_forcing_write_restart_times(forcingGroupHead)

  end subroutine ocn_time_varying_forcing_write_restart_times!}}}

!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_land_ice_forcing
!
!> \brief
!> \author Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_land_ice_forcing(domain)!{{{

    type (domain_type) :: domain

    ! create the land-ice forcing group
    call mpas_forcing_init_group(forcingGroupHead, "ocn_land_ice_forcing", domain, config_time_varying_land_ice_forcing_start_time, &
         config_time_varying_land_ice_forcing_cycle_start, config_time_varying_land_ice_forcing_cycle_duration, config_do_restart)

    call mpas_forcing_init_field(domain % streamManager, forcingGroupHead, "ocn_land_ice_forcing", &
         "landIcePressure", "land_ice_forcing", "timeVaryingForcing", "landIcePressureForcing", "linear", &
         config_time_varying_land_ice_forcing_reference_time, config_time_varying_land_ice_forcing_interval)

    call mpas_forcing_init_field(domain % streamManager, forcingGroupHead, "ocn_land_ice_forcing", &
         "landIceFraction", "land_ice_forcing", "timeVaryingForcing", "landIceFractionForcing", "linear", &
         config_time_varying_land_ice_forcing_reference_time, config_time_varying_land_ice_forcing_interval)

    call mpas_forcing_init_field(domain % streamManager, forcingGroupHead, "ocn_land_ice_forcing", &
         "landIceDraft", "land_ice_forcing", "timeVaryingForcing", "landIceDraftForcing", "linear", &
         config_time_varying_land_ice_forcing_reference_time, config_time_varying_land_ice_forcing_interval)

    call mpas_forcing_init_field_data(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_land_ice_forcing", &               ! forcingGroupName
         domain % streamManager, &               ! streamManager
         config_do_restart, &                    ! restart
         .true.)                                 ! interpolateAtInit

  end subroutine init_land_ice_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  land_ice_forcing
!
!> \brief Compute land-ice pressure
!> \author Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details Interpolate the pressure from overlying land ice
!>
!
!-----------------------------------------------------------------------

  subroutine land_ice_forcing(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    type (MPAS_time_type) :: currentForcingTime

    type (block_type), pointer :: block

!    character(len=StrKIND), pointer :: config_dt

    type (MPAS_timeInterval_type) :: timeStepESMF

    character(len=StrKIND) :: timeStamp

    real (kind=RKIND) :: dtSim

    integer :: err

    type (mpas_pool_type), pointer :: &
         mesh, &
         forcingPool, &
         timeVaryingForcingPool

    real(kind=RKIND), dimension(:), pointer :: &
         landIceFractionForcing, &
         landIceFraction, &
         landIcePressureForcing, &
         landIcePressure, &
         landIceDraftForcing, &
         landIceDraft

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    ! convert config_dt to real
    call mpas_set_timeInterval(timeStepESMF, timeString=config_dt,ierr=err)
    call mpas_get_timeInterval(timeStepESMF, dt=dtSim)

    ! use the forcing layer to get data
    call mpas_forcing_get_forcing(forcingGroupHead, "ocn_land_ice_forcing", streamManager, dtSim)

    call mpas_forcing_get_forcing_time(forcingGroupHead, "ocn_land_ice_forcing", currentForcingTime)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)
       call MPAS_pool_get_subpool(block % structs, "timeVaryingForcing", timeVaryingForcingPool)
   
       call MPAS_pool_get_dimension(mesh, "nCells", nCells)
   
       call MPAS_pool_get_array(timeVaryingForcingPool, "landIceFractionForcing", landIceFractionForcing)
       call MPAS_pool_get_array(timeVaryingForcingPool, "landIcePressureForcing", landIcePressureForcing)
       call MPAS_pool_get_array(timeVaryingForcingPool, "landIceDraftForcing", landIceDraftForcing)
       call MPAS_pool_get_array(forcingPool, "landIceFraction", landIceFraction)
       call MPAS_pool_get_array(forcingPool, "landIcePressure", landIcePressure)
       call MPAS_pool_get_array(forcingPool, "landIceDraft", landIceDraft)
   
       do iCell = 1, nCells
   
          landIceFraction(iCell) = landIceFractionForcing(iCell)
          landIcePressure(iCell) = landIcePressureForcing(iCell)
          landIceDraft(iCell) = landIceDraftForcing(iCell)
   
       enddo
       block => block % next
    end do

  end subroutine land_ice_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_aot_forcing
!
!> \brief
!> \author Collin Victor, LANL
!> \date 9/22/2021
!> \details This function is empty right now, will be updated with functionality in the future.
!>
!
!-----------------------------------------------------------------------

  subroutine init_aot_forcing(domain)!{{{

    type (domain_type) :: domain
    
    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: timeVaryingForcingPool
    real (kind=RKIND), dimension(:,:), pointer :: &
      velocityMeridionalObs, &
      velocityZonalObs, &
      normalVelocityObs
 ! 
 !   select case (trim(config_time_varying_aot_forcing_type))

 !     case ("WINDPRES")
 !        call init_aot_forcing_WINDPRES(domain)
 !     case default
 !   call mpas_log_write("aot forcing type unknown: "//trim(config_time_varying_aot_forcing_type), MPAS_LOG_CRIT)
 !   end select
    if (config_aot_obs_gen) return
 !     block => domain % blocklist
 !       do while (associated(block))
 !         call MPAS_pool_get_subpool(block % structs, 'timeVaryingForcingPool', timeVaryingForcingPool)
 !         call mpas_pool_get_array(timeVaryingForcingPool, 'velocityMeridionalObs', velocityMeridionalObs)  
 !         call mpas_pool_get_array(timeVaryingForcingPool, 'velocityZonalObs', velocityZonalObs)
 !         velocityMeridionalObs = velocityMeridional
 !         velocityZonalObs = velocityZonal
 !         block => block % next
 !       end do   

    write (*,*) "Started initializing AOT forcing."
    
    ! create the aot forcing group
    call MPAS_forcing_init_group(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_aot_forcing", &            ! forcingGroupName
         domain, &                               ! domain
         config_time_varying_aot_forcing_start_time, &            ! startTimeStr
         config_time_varying_aot_forcing_cycle_start, &           ! forcingCycleStart
         config_time_varying_aot_forcing_cycle_duration, &        ! forcingCycleDuration
         config_do_restart)!, &                  ! restart
         !.false.)                               ! forcingCycleStartInclusive
    
    select case (trim(config_time_varying_aot_data_type))
    case ('edges')
    ! normalVelocityObs init
    call MPAS_forcing_init_field(&
         domain % streamManager, &               ! streamManager
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_aot_forcing", &            ! forcingGroupName
         "normalVelocityObs", &                         ! forcingName
         "aot_forcing_edges", &                ! forcingStreamID
         "timeVaryingForcing", &                 ! poolname
         "normalVelocityObs", &                         ! fieldname
         config_time_varying_aot_forcing_type, &                             ! interpolationType
         config_time_varying_aot_forcing_reference_time, &                 ! forcingReferenceTimeStr
         config_time_varying_aot_forcing_interval, &                      ! forcingIntervalStr
         "next")                                ! forcingInitializationType

    case ('cells')
      ! velocityZonalObs init
      call MPAS_forcing_init_field(&
           domain % streamManager, &               ! streamManager
           forcingGroupHead, &                     ! forcingGroupHead
           "ocn_aot_forcing", &            ! forcingGroupName
           "velocityZonalObs", &                         ! forcingName
           "aot_forcing_cells", &                ! forcingStreamID
           "timeVaryingForcing", &                 ! poolname
           "velocityZonalObs", &                         ! fieldname
           "linear", &                             ! interpolationType
           config_time_varying_aot_forcing_reference_time, &                 ! forcingReferenceTimeStr
           config_time_varying_aot_forcing_interval)!,                      ! forcingIntervalStr
           !"next")                                ! forcingInitializationType
  
      ! velocityMeridionalObs init
      call MPAS_forcing_init_field(&
           domain % streamManager, &
           forcingGroupHead, &
           "ocn_aot_forcing", &
           "velocityMeridionalObs", &
           "aot_forcing_cells", &
           "timeVaryingForcing", &
           "velocityMeridionalObs", &
           "linear", &
           config_time_varying_aot_forcing_reference_time, &
           config_time_varying_aot_forcing_interval)!, &
           !"next")
    
    case ('both')
    ! velocityNormalObs init 
    call MPAS_forcing_init_field(&
         domain % streamManager, &               ! streamManager
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_aot_forcing", &            ! forcingGroupName
         "normalVelocityObs", &                         ! forcingName
         "aot_forcing_edges", &                ! forcingStreamID
         "timeVaryingForcing", &                 ! poolname
         "normalVelocityObs", &                         ! fieldname
         config_time_varying_aot_forcing_type, &                             ! interpolationType
         config_time_varying_aot_forcing_reference_time, &                 ! forcingReferenceTimeStr
         config_time_varying_aot_forcing_interval, &                      ! forcingIntervalStr
         "next")                                ! forcingInitializationType
      ! velocityZonalObs init
      call MPAS_forcing_init_field(&
           domain % streamManager, &               ! streamManager
           forcingGroupHead, &                     ! forcingGroupHead
           "ocn_aot_forcing", &            ! forcingGroupName
           "velocityZonalObs", &                         ! forcingName
           "aot_forcing_cells", &                ! forcingStreamID
           "timeVaryingForcing", &                 ! poolname
           "velocityZonalObs", &                         ! fieldname
           "linear", &                             ! interpolationType
           config_time_varying_aot_forcing_reference_time, &                 ! forcingReferenceTimeStr
           config_time_varying_aot_forcing_interval)!,                      ! forcingIntervalStr
           !"next")                                ! forcingInitializationType
  
      ! velocityMeridionalObs init
      call MPAS_forcing_init_field(&
           domain % streamManager, &
           forcingGroupHead, &
           "ocn_aot_forcing", &
           "velocityMeridionalObs", &
           "aot_forcing_cells", &
           "timeVaryingForcing", &
           "velocityMeridionalObs", &
           "linear", &
           config_time_varying_aot_forcing_reference_time, &
           config_time_varying_aot_forcing_interval)!, &
           !"next")
    
    
    end select

    
    if (config_use_aot_salt) then
    ! salinityObs init
    call MPAS_forcing_init_field(&
         domain % streamManager, &               ! streamManager
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_aot_forcing", &            ! forcingGroupName
         "salinityObs", &                         ! forcingName
         "aot_forcing_salt", &                ! forcingStreamID
         "timeVaryingForcing", &                 ! poolname
         "salinityObs", &                         ! fieldname
         config_time_varying_aot_forcing_type, &                             ! interpolationType
         config_time_varying_aot_forcing_reference_time, &                 ! forcingReferenceTimeStr
         config_time_varying_aot_forcing_interval, &                      ! forcingIntervalStr
         "next")                                ! forcingInitializationType

    endif
    if (config_use_aot_temp) then
    ! salinityObs init
    call MPAS_forcing_init_field(&
         domain % streamManager, &               ! streamManager
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_aot_forcing", &            ! forcingGroupName
         "temperatureObs", &                         ! forcingName
         "aot_forcing_temp", &                ! forcingStreamID
         "timeVaryingForcing", &                 ! poolname
         "temperatureObs", &                         ! fieldname
         config_time_varying_aot_forcing_type, &                             ! interpolationType
         config_time_varying_aot_forcing_reference_time, &                 ! forcingReferenceTimeStr
         config_time_varying_aot_forcing_interval, &                      ! forcingIntervalStr
         "next")                                ! forcingInitializationType
    endif

    call MPAS_forcing_init_field_data(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_aot_forcing", &            ! forcingGroupName
         domain % streamManager, &               ! streamManager
         config_do_restart, &                    ! restart
         .false.)                                ! interpolateAtInit
    ! aotPressure init
    !write (*,*) "Finished initializing AOT forcing."
  end subroutine init_aot_forcing!}}}




!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  aot_forcing
!
!> \brief Adds forcing from aot data.
!> \author Collin Victor, LANL
!> \date 9/22/2021
!> \details
!>
!
!-----------------------------------------------------------------------


  subroutine aot_forcing(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    type (MPAS_time_type) :: currentForcingTime

    type (block_type), pointer :: block

    type (MPAS_timeInterval_type) :: timeStepESMF

    character(len=StrKIND) :: timeStamp

    real (kind=RKIND) :: dtSim

    integer :: err
    
    !don't try to update data from stream if disabled
    if (config_aot_obs_gen) return
    ! convert config_dt to real
    call mpas_set_timeInterval(timeStepESMF, timeString=config_dt,ierr=err)
    call mpas_get_timeInterval(timeStepESMF, dt=dtSim)

    ! use the forcing layer to get data
    call MPAS_forcing_get_forcing(&
         forcingGroupHead, &               ! forcingGroupHead
         "ocn_aot_forcing", &      ! forcingGroupName
         streamManager, &                  ! streamManager
         dtSim)                            ! dt

    call MPAS_forcing_get_forcing_time(&
         forcingGroupHead, &               ! forcingGroupHead
         "ocn_aot_forcing", &      ! forcingGroupName
         currentForcingTime)               ! forcingTime
    
!    write (*,*) "Updated AOT forcing."
    !call mpas_get_time(curr_time=currentForcingTime, dateTimeString=timeStamp, ierr=err)
    !call mpas_log_write('Forcing time ' // trim(timeStamp))

 !!! Uncomment and modify this section to add any modifications t the aot forcing after loading.
   ! ! perform post forcing
   ! block => domain % blocklist
   ! do while (associated(block))
   !
   !    ! apply drag law to compute wind stress from wind speed
   !    if (trim(config_time_varying_aot_forcing_type) == "WINDPRES") then
   !       call post_aot_forcing(block)
   !    endif
   !
   !    block => block % next
   ! end do
  end subroutine aot_forcing!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_aot_obs_update
!
!> \brief Generates observational data to save to file
!> \author Collin Victor, LANL
!> \date 9/22/2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_aot_obs_update_vel(timeVaryingForcingPool, statePool)!{{{
    !------------------------------------------------------------------
    ! input variables
    !------------------------------------------------------------------
    type (mpas_pool_type), intent(in) :: timeVaryingForcingPool 
    type (mpas_pool_type), intent(in) :: statePool 

    !------------------------------------------------------------------
    ! local variables
    !------------------------------------------------------------------

    real (kind=RKIND), dimension(:,:), pointer :: &
      velocityMeridionalObs,  & !< [in] meridional velocity at cell centers
      velocityZonalObs,       & !< [in] zonal velocity at cell centers    
      normalVelocityObs,      & !< [in] normal velocity at edges
      normalVelocity            !< [in] normal velocity at edges

 
  if (.not. config_aot_obs_gen) return
  select case (trim(config_time_varying_aot_data_type))
  case ('edges')
    call mpas_pool_get_array(timeVaryingForcingPool,'normalVelocityObs', normalVelocityObs)
    call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
    normalVelocityObs = normalVelocity
  case ('cells')
    call mpas_pool_get_array(timeVaryingForcingPool,'velocityZonalObs', velocityZonalObs)
    call mpas_pool_get_array(timeVaryingForcingPool,'velocityMeridionalObs', velocityMeridionalObs)
    velocityZonalObs = velocityZonal
    velocityMeridionalObs = velocityMeridional
  case ('both')
    call mpas_pool_get_array(timeVaryingForcingPool,'normalVelocityObs', normalVelocityObs)
    call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
    normalVelocityObs = normalVelocity
    call mpas_pool_get_array(timeVaryingForcingPool,'velocityZonalObs', velocityZonalObs)
    call mpas_pool_get_array(timeVaryingForcingPool,'velocityMeridionalObs', velocityMeridionalObs)
    velocityZonalObs = velocityZonal
    velocityMeridionalObs = velocityMeridional
  end select
  end subroutine ocn_time_varying_forcing_aot_obs_update_vel!}}}



  subroutine ocn_time_varying_forcing_aot_obs_update_tracers(timeVaryingForcingPool,tracers, saltInd, tempInd, nCellsSolve, minLevelCell, maxLevelCell, nVertLevels)!{{{
    !------------------------------------------------------------------
    ! input variables
    !------------------------------------------------------------------
    real (kind=RKIND), dimension(:,:,:), pointer,  intent(in) :: tracers
    integer, dimension(:), intent(in) :: minLevelCell, maxLevelCell
    integer, intent(in) :: saltInd, tempInd, nCellsSolve, nVertLevels
    type (mpas_pool_type), intent(inout) :: timeVaryingForcingPool
  
    !------------------------------------------------------------------
    ! local variables
    !------------------------------------------------------------------
    integer :: iCell, iLevel
    real (kind=RKIND), dimension(:,:), pointer :: &
      salinityObs,  & !< [in] meridional velocity at cell centers
      temperatureObs  !< [in] zonal velocity at cell centers    

    real(kind=RKIND), target :: saltCell, tempCell

!  write(*,*) "Started AOT observation generator routine for tracers" 
  if (.not. config_aot_obs_gen) return

!  write(*,*) "Test salinity:"
!  write(*,*) saltInd 
  if (associated(tracers)) then
!    write(*,*) "Tracer is associated"  
!  call mpas_log_write("Testing Salinity: ", MPAS_LOG_WARN, realArgs = (/ tracers(saltInd, 2, 2)))
!    if (config_use_aot_salt) then
      call mpas_pool_get_array(timeVaryingForcingPool, 'salinityObs', salinityObs)
!      write(*,*) salinityObs(1,1)
      salinityObs = tracers(saltInd,:,:)
!      write(*,*) salinityObs(1,1)
  !    do iCell=1,nCellsSolve
  !      do iLevel=minLevelCell(iCell), maxLevelCell(iCell)
  !!        saltCell = tracers(saltInd,iLevel,iCell)
  !!        salinityObs(iLevel, iCell) => saltCell
  !        salinityObs(iLevel, iCell) = tracers(saltInd,iLevel,iCell)
  !      enddo
  !    enddo
    endif
!    write(*,*) "Test temperature:"
!    write(*,*) tracers(tempInd, 2,2)
    if (config_use_aot_temp) then
      call mpas_pool_get_array(timeVaryingForcingPool, 'temperatureObs', temperatureObs)
  !        tempCell = tracers(tempInd,iLevel,iCell)
  !        temperatureObs(iLevel, iCell) => tempCell
          temperatureObs = tracers(tempInd,:,:)
    endif
!  endif
  end subroutine ocn_time_varying_forcing_aot_obs_update_tracers!}}}



end module ocn_time_varying_forcing
! vim: foldmethod=marker
