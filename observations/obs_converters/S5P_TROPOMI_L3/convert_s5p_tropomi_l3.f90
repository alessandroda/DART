program convert_s5p_tropomi_L2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! author: Alessandro D'Ausilio
! email: a.dausilio@aria-net.it
! convert_S5p_tropomi_L2 - reads the satellite data as downloaded from the
! CAMS satellite operator.
! It sets also the metadata to be included in the obs_seq.out in order to execute
! the forward operator.
! 07.03.2024 -> add metadata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use types_mod, only : r8, missing_r8
   use location_mod, only : VERTISPRESSURE
   use netcdf
   use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
   use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
      static_init_obs_sequence, init_obs, write_obs_seq, &
      init_obs_sequence, get_num_obs, &
      set_copy_meta_data, set_qc_meta_data
   use obs_kind_mod, only : SAT_NO2_TROPOMI
   use obs_utilities_mod
   use sort_mod, only : index_sort
   use sat_obs_mod,   only : T_SatObs, ReadSatObs, SatObsDone
   use time_manager_mod, only : time_type, set_calendar_type, set_date, &
      increment_time, get_time, operator(-), GREGORIAN
   use utilities_mod, only : initialize_utilities, finalize_utilities, find_namelist_in_file, check_namelist_read, &
      nmlfileunit, do_nml_file, do_nml_term
   use obs_def_SAT_NO2_TROPOMI_mod, only : set_obs_def_tropomi

   implicit none

   type datetime
      integer :: year
      integer :: month
      integer :: day
      integer :: hour
      integer :: minute
      integer :: second
      integer :: millisecond
   end type datetime

   type(obs_sequence_type) :: obs_seq
   type(obs_type)          :: obs, prev_obs
   type(T_SatObs)          :: tsat_obs
   type(time_type)         :: comp_day0, time_obs, prev_time, gregorian_sat_obs_time
   character(len=16),  parameter ::s5p_netcdf_file = 'S5p_NO2_16742.nc'
   character(len=129), parameter :: s5p_out_file    = 'obs_seq.out'
   character(len=*),parameter        :: pollutant='NO2'
   integer,  allocatable :: used(:), tused(:), sorted_used(:)
   real(r8) :: qc, obsv, vval
   real(r8), allocatable :: lat(:), lon(:), pres(:, :), qa_value(:), amf(:), tmp(:), vcd(:), vcd_errvar(:), time(:), tobs(:)
   integer, parameter :: num_qc  = 1, num_copies = 1
   integer  :: ncid, nobs, n, i, oday, osec, nused, nlayeri, nlayer, obsindx
   integer  :: iunit, rcio ! integers to read namelist
   logical  :: file_exist, first_obs

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

   real :: dom_west = -180.0
   real :: dom_east = 180.0
   real :: dom_north = 90.0
   real :: dom_south = -90.0
   integer  :: nz = 16 ! FARM model layers
   real     :: dlon = 0.1
   real     :: dlat = 0.1
   logical  :: filter_on_model_grid = .false.
   namelist /model_grid_nml/ dom_west, dom_east, dom_south, dom_north, nz, dlon, dlat, filter_on_model_grid

   call initialize_utilities('convert_s5p_tropomi_l3')


! Read the namelist entry
   call find_namelist_in_file("input.nml", "model_grid_nml", iunit)
   read(iunit, nml = model_grid_nml, iostat = rcio)
   call check_namelist_read(iunit, rcio, "model_grid_nml")

   ! Record the namelist values used for the run ...
   if (do_nml_file()) write(nmlfileunit, nml=model_grid_nml)
   if (do_nml_term()) write(     *     , nml=model_grid_nml)
   ! put the reference date into DART format
   call set_calendar_type(GREGORIAN)
   first_obs = .true.

   !  either read existing obs_seq or create a new one
   call static_init_obs_sequence()
   call init_obs(obs,      num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)

   inquire(file=s5p_out_file, exist=file_exist)

   if ( file_exist ) then

      ! existing file found, append to it
      call read_obs_seq(s5p_out_file, 0, 0, 2*nobs, obs_seq)

   else
      call tsat_obs%Read(s5p_netcdf_file, pollutant)
      nobs = SIZE(tsat_obs%vcd)
      allocate(used(nobs))
      allocate(tused(nobs))
      allocate(tobs(nobs))
      allocate(sorted_used(nobs))
      ! create a new one
      call init_obs_sequence(obs_seq, num_copies, num_qc, 2*nobs)
      do i = 1, num_copies
         call set_copy_meta_data(obs_seq, i, 'S5P observation')
      end do
      do i = 1, num_qc
         call set_qc_meta_data(obs_seq, i, 'Data QC')
      end do

   endif

   ! Set the DART data quality control.  Be consistent with NCEP codes;
   ! 0 is 'must use', 1 is good, no reason not to use it.
   qc = 1.0_r8
   nused = 0
   obsloop1: do n = 1, nobs

      ! check the lat/lon values to see if they are ok
      if ( tsat_obs%lat(n) >  90.0_r8 .or. tsat_obs%lat(n) <  -90.0_r8 ) cycle obsloop1
      if ( tsat_obs%lon(n) > 180.0_r8 .or. tsat_obs%lon(n) < -180.0_r8 ) cycle obsloop1
      if (filter_on_model_grid) then
         ! Check if lat/lon values are within the specified domain
         if ( tsat_obs%lat(n) > dom_north .or. tsat_obs%lat(n) < dom_south ) cycle obsloop1
         if ( tsat_obs%lon(n) > dom_east .or. tsat_obs%lon(n) < dom_west ) cycle obsloop1
      endif
      ! change lon from -180 to 180 into 0-360
      if ( tsat_obs%lon(n) < 0.0_r8 )  tsat_obs%lon(n) = tsat_obs%lon(n) + 360.0_r8
      ! filter obs with qa lower than 0.75
      if (tsat_obs%qa_flag(n) < 0.75) cycle obsloop1
      ! filter sat observations negative
      if (tsat_obs%vcd(1, n) < 0) cycle obsloop1
      ! the 'used' array are the index numbers of used obs
      ! the 'tused' array are the times of those obs so we can
      ! sort them later by time.
      nused = nused + 1
      used(nused) = n
   end do obsloop1

   ! sort by time
   call index_sort(tused, sorted_used, nused)

   obsloop2: do i = 1, nused

      ! get the next unique observation in sorted time order
      n = used(sorted_used(i))

      ! compute time of observation
      time_obs = set_date(tsat_obs%date_time(i)%year,tsat_obs%date_time(i)%month,tsat_obs%date_time(i)%day, &
         tsat_obs%date_time(i)%hour, tsat_obs%date_time(i)%minute, tsat_obs%date_time(i)%second)

      ! extract actual time of observation in file into oday, osec.
      call get_time(time_obs, osec, oday)

      call set_obs_def_tropomi(n, REAL(tsat_obs%kernel_trop(1,:,n), 8), REAL(tsat_obs%pressure(:,n), 8))

      call create_3d_obs(REAL(tsat_obs%lat(n),8),REAL(tsat_obs%lon(n), 8), 1.0_r8, VERTISPRESSURE, REAL(tsat_obs%vcd(1, n), 8), &
         SAT_NO2_TROPOMI, REAL(tsat_obs%vcd_errvar(1, 1, n), 8), oday, osec, qc, obs, key = n)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   end do obsloop2

! if we added any obs to the sequence, write it now.
   if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, s5p_out_file)

   call tsat_obs%Done
! end of main program
   call finalize_utilities()
end program
