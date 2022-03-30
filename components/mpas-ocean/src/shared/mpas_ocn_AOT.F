! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_AOT
!
!> \brief MPAS ocean explicit bottom drag
!> \author Mark Petersen
!> \date   August 2017
!> \details
!>  This module contains the routine for computing
!>  tendencies from explicit bottom drag.
!
!-----------------------------------------------------------------------

module ocn_aot

   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_forcing
   use ocn_time_varying_forcing
   use ocn_diagnostics_variables
   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_aot_vel_tend, &
             ocn_aot_tracers_tend, &
             ocn_vel_aot_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      aot_vel_off, &   ! on/off switch for explicit bottom drag
      aot_salt_off, &   ! on/off switch for explicit bottom drag
      aot_temp_off      ! on/off switch for explicit bottom drag

!   character (len=StrKIND), pointer :: config_time_varying_aot_data_type
   real (kind=RKIND) :: &
      aot_coeff_vel, &          ! drag coefficient
      aot_coeff_salt, &          ! drag coefficient
      aot_coeff_temp          ! drag coefficient

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_AOT_tend
!
!> \brief   Computes tendency term from explicit bottom drag
!> \author  Mark Petersen
!> \date    15 August 2017
!> \details
!>  This routine computes the explicit bottom drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_aot_vel_tend(timeVaryingForcingPool,normalVelocity,statePool,forcingPool, tend, err) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
!      type (domain_type) intent(inout) :: domain
      type (mpas_pool_type), intent(in) :: timeVaryingForcingPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: forcingPool
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity      !< [in] normal velocity
!!         KECell,        &!< [in] kinetic energy at cell
!!         layerThickEdge  !< [in] layer thickness at edge
!      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: &
!          velocityMeridional, &!< [in] meridional velocity at cell centers
!          velocityZonal,       &!< [in] zonal velocity at cell centers 
!          velocityMeridionalObs, &!< [in] meridional velocity at cell centers
!          velocityZonalObs        !< [in] zonal velocity at cell centers 
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge,      &! loop index for edge loop
         k,          &! vertical index of lowest active layer at edge
         j,          &! 
         iSmooth,    &
         iCell,      &
         coc,        &
         eoc,        &
         nCells,     &
         cell1, cell2 ! neighbor cell addresses across edge
      real (kind=RKIND) :: aotCounterOffset
      integer, dimension(:), allocatable :: aotLayer, counter_cells, counter_layers ! tracks cell used for aot smoothing term
!      real (kind=RKIND), dimension(:,:), allocatable :: velocityZonalAOTNew, velocityMeridionalAOTNew
      real (kind=RKIND), dimension(:,:), pointer :: &
         velocityMeridionalObs, &!< [in] meridional velocity at cell centers
         velocityZonalObs,       &!< [in] zonal velocity at cell centers 
         normalVelocityObs,       &!< [in] normal velocity at edges 
         velocityMeridionalAOT, &!< [in] interpolated meridional velocity for AOT
         velocityZonalAOT        !< [in] interpolated zonal velocity for AOT 
      integer, dimension(:), pointer :: aotForcingMask
      real (kind=RKIND), dimension(:,:), allocatable :: &
!         velocityMeridionalAOT, &!< [in] interpolated meridional velocity for AOT
!         velocityZonalAOT,      &!< [in] interpolated zonal velocity for AOT 
         velocityMeridionalAOTNew, &!< [in] interpolated meridional velocity for AOT
         velocityZonalAOTNew        !< [in] interpolated zonal velocity for AOT 


      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and return if not turned on
      !*** Otherwise start timer

      err = 0
      if (aot_vel_off) return
!      if (obsGen) then
!        if (obsEdges) then 
!          call mpas_pool_get_array(timeVaryingForcingPool, 'normalVelocityObs', normalVelocityObs)
!          normalVelocityObs = normalVelocity
!        else
!          call mpas_pool_get_array(timeVaryingForcingPool, 'velocityZonalObs', velocityZonalObs)
!          call mpas_pool_get_array(timeVaryingForcingPool, 'velocityMeridionalObs', velocityMeridionalObs)
!          velocityZonalObs = velocityZonal
!          velocityMeridionalObs = velocityMeridional
!        endif
!      return 
!      endif

      !generate observational data 
      call ocn_time_varying_forcing_aot_obs_update_vel(timeVaryingForcingPool, statePool)
      
      if (config_aot_obs_gen) return


      call mpas_timer_start('AOT u-v term')
      
      ! Explicit bottom drag term:
      ! du/dt = ... - c |u| u / h
      ! appied to bottom layer only.
      ! This term comes from the bottom boundary condition in the vertical
      ! momentum mixing, and is explicit if both |u| and u are chosen to be at
      ! time level n.
           

      select case (trim(config_time_varying_aot_data_type))
      case ('edges')
        call mpas_pool_get_array(timeVaryingForcingPool, 'normalVelocityObs', normalVelocityObs)
  
  
!      call mpas_log_write("Testing observational data loading: $r", MPAS_LOG_WARN, realArgs=(/ velocityZonalObs(1,1,1)))
#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(cellsOnEdge, tend, latCell,&
        !$acc            velocityMeridionalAOT, velocityZonalAOT) &
        !$acc    private(k)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(,)
#endif
        do iEdge = 1, nEdgesOwned
!          cell1 = cellsOnEdge(1,iEdge)
!          cell2 = cellsOnEdge(2,iEdge)
  !        k =  maxLevelEdgeTop(iEdge)
!          if (latCell(cell1)<=1.2217 .AND. latCell(cell2) <=1.2217) then
            do k = 1, maxLevelEdgeTop(iEdge)
  !        if (k > 0) then
  !           tend(k,iEdge) = tend(k,iEdge) - aot_coeff_vel* &
  !                           sqrt(KECell(k,cell1) + KECell(k,cell2))* &
  !                           normVelocity(k,iEdge)/layerThickEdge(k,iEdge)
  !        end if
              tend(k,iEdge) = tend(k,iEdge) + aot_coeff_vel*(normalVelocityObs(k,iEdge) -  normalVelocity(k,iEdge))
            enddo
        enddo
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif
  
   !--------------------------------------------------------------------

      case('cells') 
  
        call mpas_pool_get_array(timeVaryingForcingPool, 'velocityZonalObs', velocityZonalObs)
        call mpas_pool_get_array(timeVaryingForcingPool, 'velocityMeridionalObs', velocityMeridionalObs)
        
        call mpas_pool_get_array(forcingPool, 'velocityZonalAOT', velocityZonalAOT)
        call mpas_pool_get_array(forcingPool, 'velocityMeridionalAOT', velocityMeridionalAOT)
        
        call mpas_pool_get_array(forcingPool, 'aotForcingMask', aotForcingMask)
  
  
  !      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        nCells = nCellsAll
  !      allocate(velocityZonalAOTNew(nVertLevels,nCells))
  !      allocate(velocityMeridionalAOTNew(nVertLevels,nCells))
        allocate(counter_layers(nVertLevels))
        allocate(counter_cells(nCells))
        allocate(aotLayer(nCells))
  !      allocate(velocityMeridionalAOT(nVertLevels,nCells))
  !      allocate(velocityZonalAOT(nVertLevels,nCells))
        allocate(velocityMeridionalAOTNew(nVertLevels,nCells))
        allocate(velocityZonalAOTNew(nVertLevels,nCells))
  
  
  !initialize aotLayer for non-masked cells
        do iCell = 1, nCells
             counter_cells(iCell) = 0
          if(aotForcingMask(iCell) == 1) then
             counter_cells(iCell) = 0
             aotLayer(iCell) = 1 !non-masked cells are in the first layer, which increments by 1 going outward to neighbors
             do k = 1, maxLevelCell(iCell)
             velocityMeridionalAOT(k,iCell) = velocityMeridionalObs(k,iCell) - velocityMeridional(k,iCell)
             velocityZonalAOT(k,iCell) = velocityZonalObs(k,iCell) - velocityZonal(k,iCell)
                !velocityZonalAOT(k,iCell) = 1
                !velocityMeridionalAOT(k,iCell) = iCell
             enddo !k
             do k = maxLevelCell(iCell) + 1, nVertLevels
               velocityMeridionalAOT(k,iCell) = 0
               velocityZonalAOT(k,iCell) = 0
             enddo
          else
            aotLayer(iCell) = 0
              do k = 1, nVertLevels
               velocityMeridionalAOT(k,iCell) = 0
               velocityZonalAOT(k,iCell) = 0
             enddo
          endif !aotForcingMask(iCell)
        enddo !iCell
  
  
        ! Nearest Neighbor initializeation
        if(config_aot_interp_weight == "constant") then
          aotCounterOffset = 1
          do iSmooth = 1, config_aot_neighbor_iterations
            do iCell = 1, nCells
              if (aotLayer(iCell) == iSmooth) then
                do j = 1, nEdgesOnCell(iCell)
                  coc = cellsOnCell(j,iCell)
                  
                  if (coc < nCells+1) then
                    if (aotLayer(coc) == 0 .OR. aotLayer(coc) == iSmooth + 1) then
                      aotLayer(coc) = iSmooth + 1
                      !counter_cells = counter_cells + 1
                      do k = 1, min(maxLevelCell(iCell),maxLevelCell(coc)) 
                        velocityMeridionalAOT(k,coc) = (velocityMeridionalAOT(k,coc)*counter_cells(coc) &
                                                     + velocityMeridionalAOT(k,iCell))/(counter_cells(coc)+aotCounterOffset)
                        velocityZonalAOT(k,coc) = (velocityZonalAOT(k,coc)*counter_cells(coc) &
                                                + velocityZonalAOT(k,iCell))/(counter_cells(coc)+aotCounterOffset)
                      enddo !k
                      counter_cells(coc) = counter_cells(coc) + aotCounterOffset
                    endif !aotLayer(coc)
                  endif !coc
                enddo !j
              endif !aotLayer(iCell)
            enddo !iCell
            ! HALO update here if more than 3 iterations
          enddo !iSmooth
        else if(config_aot_interp_weight == "geometric") then
          do iSmooth = 1, config_aot_neighbor_iterations
            do iCell = 1, nCells
              if (aotLayer(iCell) == iSmooth) then
                do j = 1, nEdgesOnCell(iCell)
                  coc = cellsOnCell(j,iCell)
                  eoc = edgesOnCell(j,iCell)
                  if (coc < nCells+1) then
                    if (aotLayer(coc) == 0 .OR. aotLayer(coc) == iSmooth + 1) then
                      aotLayer(coc) = iSmooth + 1
                      aotCounterOffset = dcEdge(eoc)
                      !counter_cells = counter_cells + 1
                      do k = 1, maxLevelCell(iCell)
                        velocityMeridionalAOT(k,coc) = (velocityMeridionalAOT(k,coc)*counter_cells(coc) &
                                                     + velocityMeridionalAOT(k,iCell))/(counter_cells(coc)+aotCounterOffset)
                        velocityZonalAOT(k,coc) = (velocityZonalAOT(k,coc)*counter_cells(coc) &
                                                + velocityZonalAOT(k,iCell))/(counter_cells(coc)+aotCounterOffset)
                      enddo !k
                      counter_cells(coc) = counter_cells(coc) + aotCounterOffset
                    endif !aotLayer(coc)
                  endif !coc
                enddo !j
              endif !aotLayer(iCell)
            enddo !iCell
          enddo !iSmooth
        endif
  
  
        velocityZonalAOTNew  = velocityZonalAOT
        velocityMeridionalAOTNew = velocityMeridionalAOT
  ! Initialize vertical counter
  
        ! Once data has been initialized apply iterated averaging to smooth it
        if(config_aot_interp_weight == "constant") then
          aotCounterOffset = 1
          do iSmooth = 1, config_aot_smoothing_iterations
            do iCell = 1, nCells
              if (aotForcingMask(iCell) == 0) then
                do k = 1, maxLevelCell(iCell)
                  velocityZonalAOTNew(k,iCell) = velocityZonalAOT(k,iCell)
                  velocityMeridionalAOTNew(k,iCell) = velocityMeridionalAOT(k,iCell)
                enddo !k
                counter_layers = 1
                do j = 1, nEdgesOnCell(iCell)
                  coc = cellsOnCell(j,iCell)
                  if (coc < nCells+1 ) then
                      do k = 1, min(maxLevelCell(iCell), maxLevelCell(coc))
                          velocityMeridionalAOTNew(k,iCell) = velocityMeridionalAOTNew(k,iCell) + velocityMeridionalAOT(k,coc)
                          velocityZonalAOTNew(k,iCell) = velocityZonalAOTNew(k,iCell) + velocityZonalAOT(k,coc)
                          counter_layers(k) = counter_layers(k) + aotCounterOffset
                      enddo !k
                  endif !coc
                enddo !j
                do k = 1, maxLevelCell(iCell)
                  velocityMeridionalAOTNew(k,iCell) = velocityMeridionalAOTNew(k,iCell)/counter_layers(k)
                  velocityZonalAOTNew(k,iCell) = velocityZonalAOTNew(k,iCell)/counter_layers(k)
                enddo !k
              endif !aotLayer(iCell)
            enddo !iCell
            velocityMeridionalAOT = velocityMeridionalAOTNew
            velocityZonalAOT = velocityZonalAOTNew
          enddo !iSmooth
  
  !! replace aotCounter in smoothing with aotCounter_depth that is reset 
  
        else if(config_aot_interp_weight == "geometric") then
          do iSmooth = 1, config_aot_smoothing_iterations
            do iCell = 1, nCells
              if (aotForcingMask(iCell) == 0) then
                do k = 1, maxLevelCell(iCell)
                  velocityZonalAOTNew(k,iCell) = velocityZonalAOT(k,iCell)
                  velocityMeridionalAOTNew(k,iCell) = velocityMeridionalAOT(k,iCell)
                enddo !k
                counter_layers(:) = 0
                do j = 1, nEdgesOnCell(iCell)
                  coc = cellsOnCell(j,iCell)
                  eoc = edgesOnCell(j,iCell)
                  if (coc < nCells+1) then
                    aotCounterOffset = dcEdge(eoc)
                    do k = 1, min(maxLevelCell(coc),maxLevelCell(iCell))
                      velocityMeridionalAOTNew(k,iCell) = velocityMeridionalAOTNew(k,iCell) + velocityMeridionalAOT(k,coc)
                      velocityZonalAOTNew(k,iCell) = velocityZonalAOTNew(k,iCell) + velocityZonalAOT(k,coc)
                      counter_layers(k) = counter_layers(k) + aotCounterOffset
                    enddo !k
                  endif !coc
                enddo !j
                do k = 1, maxLevelCell(iCell)
                  velocityMeridionalAOTNew(k,iCell) = velocityMeridionalAOTNew(k,iCell)/counter_layers(k)
                  velocityZonalAOTNew(k,iCell) = velocityZonalAOTNew(k,iCell)/counter_layers(k)
                enddo !k
              endif !aotLayer(iCell)
            enddo !iCell
            velocityMeridionalAOT = velocityMeridionalAOTNew
            velocityZonalAOT = velocityZonalAOTNew
          enddo !iSmooth
        endif
  
  
  
  
  
  !      call mpas_log_write("Testing observational data loading: $r", MPAS_LOG_WARN, realArgs=(/ velocityZonalObs(1,1,1)))
#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(cellsOnEdge, tend, latCell,&
        !$acc            velocityMeridionalAOT, velocityZonalAOT) &
        !$acc    private(k, cell1, cell2)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(, cell1, cell2)
#endif
        do iEdge = 1, nEdgesOwned
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
  !        k =  maxLevelEdgeTop(iEdge)
          if (latCell(cell1)<=1.2217 .AND. latCell(cell2) <=1.2217) then
            do k = 1, maxLevelEdgeTop(iEdge)
  !        if (k > 0) then
  !           tend(k,iEdge) = tend(k,iEdge) - aot_coeff_vel* &
  !                           sqrt(KECell(k,cell1) + KECell(k,cell2))* &
  !                           normVelocity(k,iEdge)/layerThickEdge(k,iEdge)
  !        end if
              tend(k,iEdge) = tend(k,iEdge) + aot_coeff_vel*&
                (  sin(angleEdge(iEdge)) * (0.5*(velocityMeridionalAOT(k,cell1) + velocityMeridionalAOT(k,cell2))) + &
                   cos(angleEdge(iEdge)) * (0.5*(velocityZonalAOT(k,cell1) + velocityZonalAOT(k,cell2)))     )
            enddo
          end if
        enddo
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif
  
        deallocate( aotLayer, counter_cells, counter_layers, velocityZonalAOTNew, velocityMeridionalAOTNew)
      end select
      call mpas_timer_stop('AOT u-v term')

   !--------------------------------------------------------------------

   end subroutine ocn_aot_vel_tend!}}}

   subroutine ocn_aot_tracers_tend(timeVaryingForcingPool, activeTracers, saltInd, tempInd, nCellsSolve, minLevelCell, maxLevelCell, nVertLevels, tracerTend, err) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, dimension(:), intent(in) :: minLevelCell, maxLevelCell
      integer, intent(in) :: saltInd, tempInd, nCellsSolve, nVertLevels
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tracerTend

      real (kind=RKIND), dimension(:,:,:), pointer, intent(in) ::  activeTracers
      type (mpas_pool_type), intent(inout) :: timeVaryingForcingPool
      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: Error flag 
 
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, iLevel

      real (kind=RKIND), dimension(:,:), pointer :: &
         salinityObs,         &!< [in] salinity at cell centers 
         temperatureObs        !< [in] temperature at cell centers 
!      write(*,*) "Started Tracers AOT Routine" 
      err = 0
      if (aot_vel_off) return
      !generate observational data 
      call ocn_time_varying_forcing_aot_obs_update_tracers(timeVaryingForcingPool, activeTracers, saltInd, tempInd, nCellsSolve, minLevelCell, maxLevelCell, nVertLevels)
      
      if (config_aot_obs_gen) return
!      write(*,*) "Starting AOT term calculation"
      if (config_use_aot_salt) then
        call mpas_pool_get_array(timeVaryingForcingPool, 'salinityObs', salinityObs)
        do iCell=1,nCellsSolve
          do iLevel=minLevelCell(iCell), maxLevelCell(iCell)
            tracerTend(saltInd, iLevel, iCell) = tracerTend(saltInd,iLevel,iCell) + aot_coeff_salt*(salinityObs(iLevel,iCell) - activeTracers(saltInd, iLevel, iCell)) 
          enddo
        enddo
      endif 
      if (config_use_aot_temp) then
        call mpas_pool_get_array(timeVaryingForcingPool, 'temperatureObs', temperatureObs)
        do iCell=1,nCellsSolve
          do iLevel=minLevelCell(iCell), maxLevelCell(iCell)
            tracerTend(tempInd, iLevel, iCell) = tracerTend(tempInd,iLevel,iCell) + aot_coeff_temp*(temperatureObs(iLevel,iCell) - activeTracers(tempInd, iLevel, iCell)) 
          enddo
        enddo
      endif 
!      write(*,*) "Ended AOT term calculation"
!      write(*,*) "Ended Tracers AOT Routine"
   end subroutine ocn_aot_tracers_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_AOT_init
!
!> \brief   Initializes ocean explicit bottom drag forcing
!> \author  Mark Petersen
!> \date    August 2017
!> \details
!>  This routine initializes quantities related to explicit bottom drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_aot_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize return error code and set module defaults
      err = 0

      aot_vel_off  = .true.
      aot_salt_off  = .true.
      aot_temp_off  = .true.
      aot_coeff_vel = 0.0_RKIND
      aot_coeff_salt = 0.0_RKIND
      aot_coeff_temp = 0.0_RKIND

      !*** Reset values based on input model configuration

      if (config_use_aot_vel) then
         aot_vel_off = .false.
         aot_coeff_vel = config_aot_coeff_vel
!         write(*,*) "here"
      endif
      if (config_use_aot_salt) then
         aot_salt_off = .false.
         aot_coeff_salt = config_aot_coeff_salt
!         write(*,*) "here"
      endif
      if (config_use_aot_temp) then
         aot_temp_off = .false.
         aot_coeff_temp = config_aot_coeff_temp
!         write(*,*) "here"
      endif
!      config_time_varying_aot_data_type = config_time_varying_aot_data_type
      if (config_disable_vel_aot) then
         aot_vel_off = .true.
         aot_coeff_vel = 0.0_RKIND
!         write(*,*) "no, here"
      endif
      if (config_disable_tracers_aot) then
         aot_salt_off = .true.
         aot_temp_off = .true.
         aot_coeff_salt = 0.0_RKIND
         aot_coeff_temp = 0.0_RKIND
!         write(*,*) "no, here"
      endif
      
!      write(*,*) "Mu:"
!      write(*,*) aot_coeff_vel
   !--------------------------------------------------------------------

   end subroutine ocn_vel_aot_init!}}}

!***********************************************************************

end module ocn_aot

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
