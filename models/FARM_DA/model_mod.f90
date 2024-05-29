! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

   use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength

   use time_manager_mod, only : time_type, set_time, set_calendar_type

   use     location_mod, only : location_type, get_close_type, &
      loc_get_close_obs => get_close_obs, &
      loc_get_close_state => get_close_state, &
      set_location, set_location_missing

   use    utilities_mod, only : register_module, error_handler, &
      E_ERR, E_MSG, &
      nmlfileunit, do_output, do_nml_file, do_nml_term,  &
      find_namelist_in_file, check_namelist_read

   use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
      nc_add_global_creation_time, &
      nc_begin_define_mode, nc_end_define_mode, &
      nc_open_file_readonly, nc_get_dimension_size, &
      nc_variable_exists, nc_get_variable, nc_get_variable_size

   use state_structure_mod, only : add_domain, get_domain_size

   use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
   use default_model_mod, only : pert_model_copies, read_model_time, write_model_time, &
      init_time => fail_init_time, &
      init_conditions => fail_init_conditions, &
      convert_vertical_obs, convert_vertical_state, adv_1step

   implicit none
   private

! routines required by DART code - will be called from filter and other
! DART executables.
   public :: get_model_size,         &
      get_state_meta_data,    &
      model_interpolate,      &
      end_model,              &
      static_init_model,      &
      nc_write_model_atts,    &
      get_close_obs,          &
      get_close_state,        &
      pert_model_copies,      &
      convert_vertical_obs,   &
      convert_vertical_state, &
      read_model_time,        &
      adv_1step,              &
      init_time,              &
      init_conditions,        &
      shortest_time_between_assimilations, &
      write_model_time


   character(len=256), parameter :: source   = "model_mod.f90"
   logical :: module_initialized = .false.
   integer :: dom_id ! used to access the state structure
   type(time_type) :: assimilation_time_step

! Example Namelist
! Use the namelist for options to be set at runtime.
   character(len=256) :: farm_template_filename = './farminput.nc'
   integer  :: time_step_days      = 0
   integer  :: time_step_seconds   = 3600
   integer, parameter :: MAX_NUM_STATE_VARIABLES = 30
   integer, parameter :: MAX_NUM_COLUMNS = 6
   character(len=300) :: state_variables(MAX_NUM_STATE_VARIABLES * MAX_NUM_COLUMNS) = ' '
   namelist /model_nml/ farm_template_filename, time_step_days, time_step_seconds, state_variables

!> Metadata from the template netCDF file that describes
!> where the variable data is located and what size it is.
   type farm_1d_array
      integer  :: nsize
      real(r8), allocatable :: vals(:)
   end type


   type farm_grid
      type(farm_1d_array) :: lons
      type(farm_1d_array) :: lats
      type(farm_1d_array) :: levs
      !corner points
      real :: west
      real :: south
      real :: east
      real :: north
      !grid size
      integer :: nlon
      integer :: nlat
      integer :: nlev
      !Lon Lat grid specs
      real :: delta_lon
      real :: delta_lat
   end type farm_grid

   type(farm_grid) :: grid_data

   integer, parameter :: MAX_STATE_VARIABLES = 1
   integer, parameter :: num_state_table_columns = 5

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

   subroutine static_init_model()

      integer  :: iunit, io

      module_initialized = .true.

! Print module information to log file and stdout.
      call register_module(source)

      call find_namelist_in_file("input.nml", "model_nml", iunit)
      read(iunit, nml = model_nml, iostat = io)
      call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run
      if (do_nml_file()) write(nmlfileunit, nml=model_nml)
      if (do_nml_term()) write(     *     , nml=model_nml)

      call read_FARM_definition(farm_template_filename)

! TODO: handle time

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime
! feel free to hardcode it and remove from the namelist.
      assimilation_time_step = set_time(time_step_seconds, &
         time_step_days)


! Define which variables are in the model state
      dom_id = add_domain(farm_template_filename, num_vars=1, var_names=(/'c_NO2'/))

   end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer.

   function get_model_size()

      integer(i8) :: get_model_size

      if ( .not. module_initialized ) call static_init_model

      get_model_size = get_domain_size(dom_id)

   end function get_model_size


!------------------------------------------------------------------
! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.
!
! For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

   subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


      type(ensemble_type), intent(in) :: state_handle
      integer,             intent(in) :: ens_size
      type(location_type), intent(in) :: location
      integer,             intent(in) :: qty
      real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
      integer,            intent(out) :: istatus(ens_size)

      if ( .not. module_initialized ) call static_init_model

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
      expected_obs(:) = MISSING_R8

! istatus for successful return should be 0.
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
      istatus(:) = 1

   end subroutine model_interpolate



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

   function shortest_time_between_assimilations()

      type(time_type) :: shortest_time_between_assimilations

      if ( .not. module_initialized ) call static_init_model

      shortest_time_between_assimilations = assimilation_time_step

   end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

   subroutine get_state_meta_data(index_in, location, qty)

      integer(i8),         intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer,             intent(out), optional :: qty


      if ( .not. module_initialized ) call static_init_model

! should be set to the actual location using set_location()
      location = set_location_missing()

! should be set to the physical quantity, e.g. QTY_TEMPERATURE
      if (present(qty)) qty = 0

   end subroutine get_state_meta_data


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
   subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
      num_close, close_ind, dist, ens_handle)

      type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
      integer,                       intent(in)    :: base_type     ! observation TYPE
      type(location_type),           intent(inout) :: base_loc      ! location of interest
      type(location_type),           intent(inout) :: locs(:)       ! obs locations
      integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
      integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
      integer,                       intent(out)   :: num_close     ! how many are close
      integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
      real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
      type(ensemble_type), optional, intent(in)    :: ens_handle

      character(len=*), parameter :: routine = 'get_close_obs'

      call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
         num_close, close_ind, dist, ens_handle)

   end subroutine get_close_obs


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
   subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
      num_close, close_ind, dist, ens_handle)

      type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
      type(location_type),           intent(inout) :: base_loc     ! location of interest
      integer,                       intent(in)    :: base_type    ! observation TYPE
      type(location_type),           intent(inout) :: locs(:)      ! state locations
      integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
      integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
      integer,                       intent(out)   :: num_close    ! how many are close
      integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
      real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
      type(ensemble_type), optional, intent(in)    :: ens_handle

      character(len=*), parameter :: routine = 'get_close_state'


      call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
         num_close, close_ind, dist, ens_handle)


   end subroutine get_close_state


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

   subroutine end_model()


   end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

   subroutine nc_write_model_atts(ncid, domain_id)

      integer, intent(in) :: ncid      ! netCDF file identifier
      integer, intent(in) :: domain_id

      if ( .not. module_initialized ) call static_init_model

! put file into define mode.

      call nc_begin_define_mode(ncid)

      call nc_add_global_creation_time(ncid)

      call nc_add_global_attribute(ncid, "model_source", source )
      call nc_add_global_attribute(ncid, "model", "template")

      call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
      call nc_synchronize_file(ncid)

   end subroutine nc_write_model_atts

   subroutine fill_farm_1d_array(ncid, varname, grid_array)
      integer,            intent(in)    :: ncid
      character(len=*),   intent(in)    :: varname
      type(farm_1d_array), intent(inout) :: grid_array

      character(len=*), parameter :: routine = 'fill_farm_1d_array'

      if(nc_variable_exists(ncid, varname)) then

         call nc_get_variable_size(ncid, varname, grid_array%nsize)
         allocate(grid_array%vals(grid_array%nsize))

         call nc_get_variable(ncid, varname, grid_array%vals, routine)
      endif

   end subroutine fill_farm_1d_array

   subroutine read_FARM_definition(file_name)
      !A.D'A: based on templates it reads the netcdf for the definition
      !stores grid information for the moment very similar to SATDA_OFFLINE
      character(len=*), intent(in) :: file_name
      integer  :: ncid, DimID, TimeDimID

      character(len=*), parameter :: routine = 'read_FARM_definition'

      call error_handler(E_MSG,routine,'reading restart ['//trim(file_name)//']')

      ncid = nc_open_file_readonly(file_name, routine)

      ! longitude - FARM uses values +/- 180, DART uses values [0,360]
      grid_data%nlon = nc_get_dimension_size(ncid, 'x', routine)
      call fill_farm_1d_array(ncid, 'x', grid_data%lons)
      where (grid_data%lons%vals < 0.0_r8) grid_data%lons%vals = grid_data%lons%vals + 360.0_r8

      ! latitiude
      grid_data%nlat = nc_get_dimension_size(ncid, 'y', routine)
      call fill_farm_1d_array(ncid, 'y', grid_data%lats)

      ! level [m]
      grid_data%nlev = nc_get_dimension_size(ncid, 'z', routine)
      call fill_farm_1d_array(ncid, 'z', grid_data%levs)

      ! Get lon and lat grid specs
      grid_data%west = grid_data%lons%vals(1)
      grid_data%east= grid_data%lons%vals(grid_data%nlon)
      grid_data%south = grid_data%lats%vals(1)
      grid_data%north = grid_data%lats%vals(grid_data%nlat)
      grid_data%delta_lon      = abs((grid_data%lons%vals(1)-grid_data%lons%vals(2)))
      grid_data%delta_lat      = abs((grid_data%lats%vals(1)-grid_data%lats%vals(2)))

      ! Print all elements of grid_data
      print *, 'nlon:', grid_data%nlon
      print *, 'nlat:', grid_data%nlat
      print *, 'nlev:', grid_data%nlev
      print *, 'west:', grid_data%west
      print *, 'south:', grid_data%south
      print *, 'east:', grid_data%east
      print *, 'north:', grid_data%north
      print *, 'delta_lon:', grid_data%delta_lon
      print *, 'delta_lat:', grid_data%delta_lat

   end subroutine read_FARM_definition
end module model_mod
