! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

   use        types_mod,        only : r8, i8, MISSING_R8, MISSING_I, vtablenamelength, i4, r4

   use        time_manager_mod, only : time_type, set_time, set_calendar_type

   use        location_mod,     only : location_type, get_close_type, get_dist, query_location, get_location, &
      loc_get_close_obs => get_close_obs, &
      loc_get_close_state => get_close_state, &
      set_location, set_location_missing

   use utilities_mod, only : register_module, error_handler, &
      E_ERR, E_MSG, &
      nmlfileunit, do_output, do_nml_file, do_nml_term,              &
      find_namelist_in_file, check_namelist_read, string_to_integer, &
      string_to_real, string_to_logical, to_upper

   use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
      set_quad_coords, finalize_quad_interp, &
      quad_lon_lat_locate, quad_lon_lat_evaluate, &
      GRID_QUAD_IRREG_SPACED_REGULAR,  &
      QUAD_LOCATED_CELL_CENTERS

   use obs_kind_mod, only: QTY_NO2, QTY_PRESSURE, QTY_SURFACE_PRESSURE, get_index_for_quantity, get_num_quantities, get_name_for_quantity
   use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
      nc_add_global_creation_time, &
      nc_begin_define_mode, nc_end_define_mode, &
      nc_open_file_readonly, nc_get_dimension_size, &
      nc_variable_exists, nc_get_variable, nc_get_variable_size
   use distributed_state_mod, only : get_state, get_state_array
   use state_structure_mod, only : add_domain, get_domain_size, get_model_variable_indices, get_dart_vector_index, &
      get_num_variables

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
   integer, allocatable :: state_kinds_list(:)

! Example Namelist
! Use the namelist for options to be set at runtime.
   character(len=256) :: farm_template_filename = './farminput.nc'
   integer  :: time_step_days      = 0
   integer  :: time_step_seconds   = 3600
   integer, parameter :: MAX_NUM_STATE_VARIABLES = 30
   integer, parameter :: MAX_NUM_COLUMNS = 5
   character(len=256) :: errstring
! state_variables defines the contents of the state vector.
! each line of this input should have the form:
!
!    netcdf_variable_name, dart_quantity, clamp_min, clamp_max, update_variable
!
! all items must be strings (even if numerical values).
! for no clamping, use the string 'NA'
! to have the assimilation change the variable use 'UPDATE', else 'NO_UPDATE'

   character(len=300) :: state_variables(MAX_NUM_STATE_VARIABLES * MAX_NUM_COLUMNS) = ' '
   namelist /model_nml/           &
      farm_template_filename,        &
      time_step_days,             &
      time_step_seconds,          &
      state_variables

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

   integer, parameter :: MAX_STATE_VARIABLES = 3
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

      call get_FARM_grid(farm_template_filename)

      call set_calendar_type('GREGORIAN')

      call set_farm_variable_info(farm_template_filename, state_variables)


      assimilation_time_step = set_time(time_step_seconds, &
         time_step_days)

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
      real(r8) :: val(ens_size)
      integer  :: which_vert, status1, varid
      integer  :: four_lons(4), four_lats(4)
      integer  :: status_array(ens_size)
      real(r8) :: lon_lat_vert(3)
      real(r8) :: quad_vals(ens_size,4)
      integer(i4) :: imem, ens_ith, index_loc_lat, index_loc_lon, index_loc_vert
      integer(i8) :: index_state
      real(r8) :: cell_lat
      logical :: isfound
      integer, dimension(2) :: cell_idxs
      integer :: ivar
      type(quad_interp_handle) :: interp_handle

      real(r8) :: ps(ens_size)
      if ( .not. module_initialized ) call static_init_model

      ivar = get_varid_from_kind(qty)
      if (ivar < 1) then
         istatus = 88
         return
      endif

      lon_lat_vert = get_location(location)
      if (qty == QTY_PRESSURE .OR. qty == QTY_NO2) then
         call find_cell(lon_lat_vert(2), lon_lat_vert(1), grid_data%lats%vals, grid_data%lons%vals, cell_idxs, isfound )
         if(isfound .neqv. .true.) then
            istatus(:) = 88
            return
         end if
         call find_index_of_value(lon_lat_vert(3), index_loc_vert, 'levs', isfound)
         if(isfound .neqv. .true.) then
            istatus(:) = 89
            return
         end if
         index_state =  get_dart_vector_index(cell_idxs(2), cell_idxs(1), index_loc_vert, dom_id, ivar)
         expected_obs = get_state(index_state,state_handle)
         istatus(:) = 0
      else if (qty == QTY_SURFACE_PRESSURE) then
         call find_cell(lon_lat_vert(2), lon_lat_vert(1), grid_data%lats%vals, grid_data%lons%vals, cell_idxs, isfound )
         if(isfound .neqv. .true.) then
            istatus(:) = 90
            return
         end if
         index_state =  get_dart_vector_index(cell_idxs(2), cell_idxs(1), index_loc_vert, dom_id, ivar)
         expected_obs = get_state(index_state,state_handle)
         istatus(:) = 0
      end if
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
!
! Parameters:
!   - index_in: Integer index into the state vector.
!   - location: Output parameter for the associated location.
!   - var_type: Optional output parameter for the physical quantity.
!------------------------------------------------------------------

   subroutine get_state_meta_data(index_in, location, var_type)

      integer(i8),         intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) :: var_type
      logical, parameter    :: debug = .false.
      character(len=*), parameter :: routine = 'get_state_meta_data'
      ! Local variables
      integer :: iloc, vloc, jloc
      integer :: myvarid, myqty
      ! Ensure the module is initialized
      if ( .not. module_initialized ) call static_init_model

      ! Get variable indices and associated quantity
      call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id = myvarid, dom_id = dom_id, kind_index=myqty)
      ! Check dimensions
      if (grid_data%lons%nsize <= iloc .or. grid_data%lats%nsize <= jloc) then
         call error_handler(E_MSG, routine, text = "Invalid index: Out of bounds for grid dimensions.")
      end if

      ! Set the location using set_location()
      location = set_location(real(grid_data%lons%vals(iloc), r8), real(grid_data%lats%vals(jloc), r8), grid_data%levs%vals(vloc), 3)
      if (present(var_type)) var_type = myqty

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


   !------------------------------------------------------------------------------------------------------------------
   !> Read the data from the various FARM grid arrays
   subroutine get_FARM_grid(file_name)
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

   end subroutine get_FARM_grid

!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!> Then calls 'add_domain()' to tell the DART code which variables to
!> read into the state vector after this code returns.
!>
!>@param variable_array  the list of variables and kinds from model_mod_nml
!>@param nfields         the number of variable/Quantity pairs specified

   subroutine set_farm_variable_info(farm_template_filename, variable_array)

      character(len=*), intent(in)  :: farm_template_filename
      character(len=*), intent(in)  :: variable_array(:)

      character(len=*), parameter :: routine = 'set_farm_variable_info:'

      integer :: i, nfields
      integer, parameter :: MAX_STRING_LEN = 128

      character(len=MAX_STRING_LEN) :: varname    ! column 1, NetCDF variable name
      character(len=MAX_STRING_LEN) :: dartstr    ! column 2, DART Quantity
      character(len=MAX_STRING_LEN) :: minvalstr  ! column 3, Clamp min val
      character(len=MAX_STRING_LEN) :: maxvalstr  ! column 4, Clamp max val
      character(len=MAX_STRING_LEN) :: updatestr  ! column 5, Update output or not

      character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
      logical  :: update_list(MAX_STATE_VARIABLES)   = .FALSE.
      integer  ::   kind_list(MAX_STATE_VARIABLES)   = MISSING_I
      real(r8) ::  clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8

      nfields = 0
      ParseVariables : do i = 1, MAX_STATE_VARIABLES

         varname   = variable_array(num_state_table_columns*i-4)
         dartstr   = variable_array(num_state_table_columns*i-3)
         minvalstr = variable_array(num_state_table_columns*i-2)
         maxvalstr = variable_array(num_state_table_columns*i-1)
         updatestr = variable_array(num_state_table_columns*i  )

         if ( varname == ' ' .and. dartstr == ' ' ) exit ParseVariables ! Found end of list.

         ! if ( varname == ' ' .or.  dartstr == ' ' ) then
         !    string1 = 'model_nml:model "state_variables" not fully specified'
         !    call error_handler(E_ERR,routine,string1,source,revision,revdate)
         ! endif

         ! Make sure DART kind is valid

         ! if( get_index_for_quantity(dartstr) < 0 ) then
         !    write(string1,'(3A)') 'there is no obs_kind "', trim(dartstr), '" in obs_kind_mod.f90'
         !    call error_handler(E_ERR,routine,string1,source,revision,revdate)
         ! endif

         call to_upper(minvalstr)
         call to_upper(maxvalstr)
         call to_upper(updatestr)

         var_names(   i) = varname
         kind_list(   i) = get_index_for_quantity(dartstr)
         clamp_vals(i,1) = string_to_real(minvalstr)
         clamp_vals(i,2) = string_to_real(maxvalstr)
         update_list( i) = string_to_logical(updatestr, 'UPDATE')
         state_kinds_list = kind_list
         nfields = nfields + 1

      enddo ParseVariables

      ! if (nfields == MAX_STATE_VARIABLES) then
      !    write(string1,'(2A)') 'WARNING: There is a possibility you need to increase ', &
      !       'MAX_STATE_VARIABLES in the global variables in model_mod.f90'

      !    write(string2,'(A,i4,A)') 'WARNING: you have specified at least ', nfields, &
      !       ' perhaps more'

      !    call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
      ! endif

      ! CAM only has a single domain (only a single grid, no nests or multiple grids)

      dom_id = add_domain(farm_template_filename, nfields, var_names, kind_list, &
         clamp_vals, update_list)

   end subroutine set_farm_variable_info

   subroutine find_cell(obs_lat, obs_lon, lat, lon, cells, cell_found)
      real(8), intent(in) :: obs_lat, obs_lon    ! Array containing latitude and longitude of the location
      real(8), intent(in) :: lat(:)         ! Array of latitude values for the grid
      real(8), intent(in) :: lon(:)         ! Array of longitude values for the grid
      integer, intent(out) :: cells(2)      ! Array containing indices of the cell
      logical, intent(out) :: cell_found    ! Flag indicating if the cell is found

      integer :: i, j
      real(8) :: lat_min, lat_max, lon_min, lon_max

      ! Initialize cell_found flag
      cell_found = .false.

      ! Check if the location is within the grid boundaries
      lat_min = minval(lat)
      lat_max = maxval(lat)
      lon_min = minval(lon)
      lon_max = maxval(lon)

      if (obs_lat < lat_min .or. obs_lat > lat_max .or. &
         obs_lon < lon_min .or. obs_lon > lon_max) then
         ! Location is outside the grid boundaries
         return
      end if

      ! Iterate over grid cells to find the location
      do i = 1, size(lat) - 1
         if (obs_lat >= lat(i) .and. obs_lat <= lat(i + 1)) then
            do j = 1, size(lon) - 1
               if (obs_lon >= lon(j) .and. obs_lon <= lon(j + 1)) then
                  ! Location is within the current cell
                  cells(1) = i
                  cells(2) = j
                  cell_found = .true.
                  return
               end if
            end do
         end if
      end do

   end subroutine find_cell


   subroutine find_index_of_value(target_value, index, component, isfound)
      real(r8), intent(in) :: target_value
      character(len=*), intent(in) :: component
      integer, intent(out) :: index
      logical,intent(out) :: isfound
      integer :: i

      ! Initialize index to 0 (not found)
      index = 0

      select case (trim(component))
       case('lats')
         ! Search for target_value in vals array
         do i = 1, size(grid_data%lats%vals)
            if (grid_data%lats%vals(i) == target_value) then
               index = i
               exit
            end if
         end do
       case('lons')
         do i = 1, size(grid_data%lons%vals)
            if (grid_data%lons%vals(i) == target_value) then
               index = i
               exit
            end if
         end do
       case('levs')
         do i = 1, size(grid_data%levs%vals)
            if (grid_data%levs%vals(i) == target_value) then
               index = i
               exit
            end if
         end do
       case default
         write(*,*) 'Invalid component specified.'
      end select

      ! Display the result
      if (index /= 0) then
         ! write(*, *) 'Index of ', target_value, ' in vals(:) is ', index
         isfound = .true.
      else
         ! write(*, *) 'Value not found in vals(:)'
         isfound = .false.
      end if
   end subroutine find_index_of_value


!--------------------------------------------------------------------
!> given a kind, return what variable number it is
!--------------------------------------------------------------------
   function get_varid_from_kind(dart_kind)

      integer, intent(in) :: dart_kind
      integer             :: get_varid_from_kind

      integer :: i

      do i = 1, get_num_variables(dom_id)
         if (dart_kind == state_kinds_list(i)) then
            get_varid_from_kind = i
            return
         endif
      end do

   end function get_varid_from_kind

end module model_mod
