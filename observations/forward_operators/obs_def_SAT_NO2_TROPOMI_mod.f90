! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! ! obs_def_SAT_NO2_TROPOMI_mod
! ! author : Alessandro D'Ausilio
! ! email : a.dausilio@aria-net.it
!****************************************************************************************************
! This section defines the forward operator for assimilating S5P Tropomi NO2 in
! FARM. The 'preprocess' function utilizes this operator to incorporate appropriate
! definitions of SAT_NO2_TROPOMI in the DEFAULT_obs_def_mod.f90 template. Subsequently, it generates
! the source files obs_def_mod.f90 and obs_kind_mod.f90, crucial for filter and other DART programs.
! The DART PREPROCESS TYPE DEFINITION excludes the keyword COMMON_CODE since the observation requires
! the forward operator.
!
! The subroutine 'get_expected_SAT_NO2_TROPOMI' is employed by the filter and performs the A*V operation.
! As of 07.03.2024, the variable 'G' appears unnecessary.
!
! The 'convert_s5p_tropomi_l3' subroutine creates obs_seq.out, generating a 3D observation and adding
! kernel and pressure vectors via set_obs_def_tropomi. It subsequently calls 'write_tropomi_no2'
! to write the observations into the observation sequence file.
!
! During filter execution, the code initially passes through 'read_tropomi_no2' to read previously set
! and written Tropomi NO2 observations, performing interpolation.
!****************************************************************************************************


! BEGIN DART PREPROCESS TYPE DEFINITIONS
! SAT_NO2_TROPOMI, QTY_NO2
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_SAT_NO2_TROPOMI_mod, only : get_expected_SAT_NO2_TROPOMI, write_tropomi_no2, read_tropomi_no2, &
!   set_obs_def_tropomi
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(SAT_NO2_TROPOMI)
!       call get_expected_SAT_NO2_TROPOMI(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(SAT_NO2_TROPOMI)
!           call read_tropomi_no2(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(SAT_NO2_TROPOMI)
!           call write_tropomi_no2(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_TROPOMI
!      case(SAT_NO2_TROPOMI)
!         call set_obs_def_tropomi(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_TROPOMI

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(SAT_NO2_TROPOMI)
!         continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_SAT_NO2_TROPOMI_mod

   use typeSizes
   use        types_mod, only : r8, MISSING_R8
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
      nmlfileunit, do_nml_file, do_nml_term, &
      check_namelist_read, find_namelist_in_file
   use     location_mod, only : location_type, set_location, get_location, &
      write_location, read_location, &
      VERTISLEVEL, VERTISPRESSURE, VERTISSURFACE, VERTISHEIGHT
   use time_manager_mod, only : time_type, read_time, write_time, &
      set_time, set_time_missing
   use  assim_model_mod, only : interpolate
   use     obs_kind_mod, only : QTY_NO2, QTY_PRESSURE
   use ensemble_manager_mod,  only : ensemble_type
   use obs_def_utilities_mod, only : track_status

   implicit none
   private

   public ::  get_expected_SAT_NO2_TROPOMI, write_tropomi_no2,read_tropomi_no2, set_obs_def_tropomi

   integer, parameter               :: max_model_levs = 16

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = &
      "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   logical, save :: module_initialized = .false.

   character(len=129) :: msgstring

   real(r8), parameter :: gravity = 9.81_r8     ! gravitational acceleration (m s^-2)
   real(r8), parameter :: density = 1000.0_r8   ! water density in kg/m^3

   integer :: max_pressure_intervals = 1000   ! increase as needed
   real(r8)   :: farm_heights(16) =(/ &
      20.,   65.,  125.,  210.,  325.,  480.,  690.,  975., 1360., &
      1880., 2580., 3525., 4805., 6290., 7790., 9290./)
! default samples the atmosphere between the surface and 200 hPa
! at the model level numbers.  if model_levels is set false,
! then the default samples at 40 heights, evenly divided in
! linear steps in pressure between the surface and top.

   logical  :: model_levels = .true.        ! if true, use model levels, ignores num_pres_int
   real(r8) :: pressure_top = 20000.0       ! top pressure in pascals
   logical  :: separate_surface_level = .true.  ! false: level 1 of 3d grid is sfc
   ! true: sfc is separate from 3d grid
   integer  :: num_pressure_intervals = 30
   integer, parameter :: levels = 34 ! hardcoded
   integer, parameter :: max_obs = 100000 ! number of intervals if model_levels is F
   integer,  dimension(max_obs)   :: tropomi_nlevels
   real,dimension(levels, max_obs) :: kernel_trop_px
   real,dimension(levels, max_obs) :: pressure_px
   character(len=6), parameter :: S5Pstring = 'FO_params'

   namelist /obs_def_SAT_NO2_TROPOMI_nml/ model_levels, pressure_top,  &
      separate_surface_level, num_pressure_intervals

contains

!------------------------------------------------------------------------------
   subroutine initialize_module()

! should be called once by the other routines in this file to be sure
! the namelist has been read and any initialization code is run.

      integer :: iunit, rc

      if (module_initialized) return

      call register_module(source, revision, revdate)
      module_initialized = .true.

! Read the namelist entry
      call find_namelist_in_file("input.nml", "obs_def_SAT_NO2_TROPOMI_nml", iunit)
      read(iunit, nml = obs_def_SAT_NO2_TROPOMI_nml, iostat = rc)
      call check_namelist_read(iunit, rc, "obs_def_SAT_NO2_TROPOMI_nml")

! Record the namelist values used for the run ...
      if (do_nml_file()) write(nmlfileunit, nml=obs_def_SAT_NO2_TROPOMI_nml)
      if (do_nml_term()) write(     *     , nml=obs_def_SAT_NO2_TROPOMI_nml)

   end subroutine initialize_module


!------------------------------------------------------------------------------
   subroutine get_expected_SAT_NO2_TROPOMI(state_handle, ens_size, location, key, val, istatus)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Author: Alessandro D'Ausilio ,  Version 0: 08/03/2024
!  Model refers to FARM
!  tropomi refers to satellite data
!
!------------------------------------------------------------------------------
      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: ens_size
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: key
      real(r8),            intent(out) :: val(ens_size)
      integer,             intent(out) :: istatus(ens_size)

      integer             :: imem
      integer             :: num_levs, lev
      real(r8)            :: p_col(ens_size, max_model_levs)
      real(r8)            :: tropomi_pres_local(ens_size, 34)
      integer             :: p_col_istatus(ens_size)
      type(location_type) :: locS
      real(r8)            :: mloc(3), mloc1(3), mloc2(3)
      if ( .not. module_initialized ) call initialize_module
      val = MISSING_R8

      mloc = get_location(location)
      if (mloc(2)>90.0_r8) then
         mloc(2)=90.0_r8
      elseif (mloc(2)<-90.0_r8) then
         mloc(2)=-90.0_r8
      endif
!     For each ensemble pass the satellite pressure
      do imem = 1, ens_size
         tropomi_pres_local(imem, :) = pressure_px(:,key)
      enddo

      istatus = 0
      p_col = MISSING_R8
      lev = 1
!     FARM pressure field at pixel position
      model_levels: do
         locS = set_location(mloc(1),mloc(2),farm_heights(lev),VERTISHEIGHT)
         call interpolate(state_handle, ens_size, locS, QTY_PRESSURE, p_col(:, lev), p_col_istatus)
         if (any(p_col_istatus /= 0)) then
            p_col(:, lev) = MISSING_R8
            num_levs = lev - 1
            exit model_levels
         endif
         lev = lev + 1
      enddo model_levels


      do imem= 1, ens_size
         val(imem) =  1
      end do
      istatus = 0
   end subroutine get_expected_SAT_NO2_TROPOMI

   subroutine read_tropomi_no2(key, ifile, fileformat)
      integer, intent(out) :: key
      integer, intent(in)  :: ifile
      character(len=32)    :: fileformat

      integer              :: tropomi_nlevels_1
      real(r8)             :: tropomi_prior_1
      real(r8)             :: tropomi_psurf_1
      integer              :: keyin

      ! Dummy implementation, replace with actual code to read data from the file
      key = 0
      tropomi_nlevels_1 = 0
      tropomi_prior_1 = 0.0_r8
      tropomi_psurf_1 = 0.0_r8
      keyin = 0

      ! Print a message indicating that this is a dummy implementation
      print *, 'Dummy implementation of read_tropomi_no2 subroutine'
   end subroutine read_tropomi_no2


   subroutine write_tropomi_no2(key, ifile, fform)
      integer, intent(in) :: key
      integer, intent(in)  :: ifile
      character(len=*),  intent(in), optional :: fform
      print *, 'Dummy implementation of write subroutine'
      ! write(ifile, *) trim(S5Pstring)
      ! write(ifile, *) key
      ! write(ifile, *) pressure_px(:,key)
      ! write(ifile, *) kernel_trop_px(:,key)

   end subroutine write_tropomi_no2

   subroutine set_obs_def_tropomi(key, kernel_trop_px_1, pressure_px_1)
      integer, intent(in) :: key
      integer :: i
      character(len=32)    :: fform
      real,dimension(34) :: kernel_trop_px_1, pressure_px_1

      do i = 1, 34
         kernel_trop_px(i,key) = kernel_trop_px_1(i)
         pressure_px(i, key) = pressure_px_1(i)
      end do

   end subroutine set_obs_def_tropomi

end module obs_def_SAT_NO2_TROPOMI_mod

! END DART PREPROCESS MODULE CODE

