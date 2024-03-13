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
   use        types_mod, only : r8, MISSING_R8, r4, MISSING_R4
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
      nmlfileunit, do_nml_file, do_nml_term, ascii_file_format, &
      check_namelist_read, find_namelist_in_file
   use     location_mod, only : location_type, set_location, get_location, &
      write_location, read_location, &
      VERTISLEVEL, VERTISPRESSURE, VERTISSURFACE, VERTISHEIGHT
   use time_manager_mod, only : time_type, read_time, write_time, &
      set_time, set_time_missing
   use  assim_model_mod, only : interpolate
   use     obs_kind_mod, only : QTY_NO2, QTY_PRESSURE, QTY_SURFACE_PRESSURE
   use ensemble_manager_mod,  only : ensemble_type
   use obs_def_utilities_mod, only : track_status

   implicit none
   private

   public ::  get_expected_SAT_NO2_TROPOMI, write_tropomi_no2,read_tropomi_no2, set_obs_def_tropomi

   integer, parameter               :: max_model_levs = 16
   integer, parameter               :: max_model_p_levs = 17

! version controlled file description for error handling, do not edit
   character(len=256), parameter :: source   = &
      "$URL$"
   character(len=32 ), parameter :: revision = "$Revision$"
   character(len=128), parameter :: revdate  = "$Date$"

   logical, save :: module_initialized = .false.
   integer  :: counts1 = 0

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
   integer, parameter :: tropomi_dim = 34 ! hardcoded
   integer, parameter :: max_obs = 1000000 ! number of intervals if model_levels is F
   integer,  dimension(max_obs)   :: tropomi_nlevels
   real(r8),dimension(tropomi_dim, max_obs) :: kernel_trop_px
   real(r8),dimension(tropomi_dim, max_obs) :: pressure_px
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
      integer             :: num_levs, nz, iz, level_ith
      real(r8),allocatable                  :: model_p(:, :)
      real(r8), allocatable    :: tropomi_pres_local(:,:)
      integer             :: p_col_istatus(ens_size)
      type(location_type) :: locS
      real(r8)            :: mloc(3), mloc1(3), mloc2(3)
!Vertical for LayerAverage
      integer                           :: indexSP_kernel,indexSP_model
      real(r8),allocatable                  :: modelP_scaled(:)
      real(r8)        :: sp(ens_size)
      integer,allocatable               :: iiv(:),jjv(:),nwv
      real(r8),allocatable                  :: dxv(:),dyv(:),wwv(:)

      if ( .not. module_initialized ) call initialize_module
      val = MISSING_R8

      allocate(tropomi_pres_local(ens_size, tropomi_dim + 1))
      tropomi_pres_local = 0.0_r8

      mloc = get_location(location)
      if (mloc(2)>90.0_r8) then
         mloc(2)=90.0_r8
      elseif (mloc(2)<-90.0_r8) then
         mloc(2)=-90.0_r8
      endif
!     For each ensemble pass the satellite pressure
      do imem = 1, ens_size
         do level_ith = 1, tropomi_dim
            tropomi_pres_local(imem, level_ith) = pressure_px(level_ith,key)
         enddo
      enddo

      allocate(model_p(ens_size, max_model_p_levs))
      istatus = 0
      nz = 1
      !     FARM pressure field at pixel position
      model_levels: do
         locS = set_location(mloc(1),mloc(2),farm_heights(nz),VERTISHEIGHT)
         call interpolate(state_handle, ens_size, locS, QTY_PRESSURE, model_p(:, nz), p_col_istatus)
         if (any(p_col_istatus /= 0)) then
            model_p(:, nz) = MISSING_R8
            num_levs = nz - 1
            exit model_levels
         endif
         nz = nz + 1
      enddo model_levels
      ! nz should be at this point 16
      ! num_levs u
!     FARM surface pressure field at pixel position
      call interpolate(state_handle, ens_size, locS, QTY_SURFACE_PRESSURE, sp(:), p_col_istatus)
      if (any(p_col_istatus /= 0)) then
         sp(:) = MISSING_R8
      endif

      !_____________________________________________________
      !PRESSURE.............................................
      ! compute pressure at faces
      do iz=1, num_levs -1
         model_p(:, iz + 1)=0.5*(model_p(:, iz + 1) +model_p(:,iz))
      enddo
      model_p(:,max_model_p_levs)=0.0 ! Pa
      !First level
      model_p(:,1)=sp(:)
      !Change Units hPa-> Pa
      model_p=model_p*100
      !......................................................
      allocate( iiv(tropomi_dim*max_model_levs))
      allocate( jjv(tropomi_dim*max_model_levs))
      allocate( wwv(tropomi_dim*max_model_levs))
      allocate( dxv(max_model_levs))
      allocate( dyv(max_model_levs))
      allocate( nwv)

      call FindIndexSurfacePressure(tropomi_dim,tropomi_pres_local(1, :), indexSP_kernel)

      call FindIndexSurfacePressure(nz+1,model_p(1, :),indexSP_model)

      do imem = 1, ens_size
         model_p(imem, :) = model_p(imem, :) / model_p(imem,indexSP_model) * tropomi_pres_local(imem, indexSP_kernel)
         model_p(imem, 1)=model_p(imem, 1)+0.1
      enddo

      call ComputeWeight(max_model_levs,tropomi_dim,model_p(1, :),tropomi_pres_local(1,:),                   & !input
         iiv(:),jjv(:),dxv(:),dyv(:),nwv,wwv(:))

      do imem= 1, ens_size
         val(imem) =  1
      end do
      istatus = 0

      deallocate( iiv )
      deallocate( jjv )
      deallocate( wwv )
      deallocate( dxv )
      deallocate( dyv )
      deallocate( nwv )
      deallocate(tropomi_pres_local)
      deallocate(model_p)

   end subroutine get_expected_SAT_NO2_TROPOMI

   subroutine read_tropomi_no2(key, ifile, fform)
      integer, intent(out) :: key
      integer :: i
      integer, intent(in)  :: ifile
      character(len=*), intent(in), optional    :: fform
      character(len=32) 		:: fileformat
      real(r8), dimension(tropomi_dim)	:: avg_kernels_1
      real(r8), dimension(tropomi_dim)	:: obs_p_1
      integer 			:: keyin

      if ( .not. module_initialized ) call initialize_module

      fileformat = "ascii"   ! supply default
      if(present(fform)) fileformat = trim(adjustl(fform))


      avg_kernels_1(:) = 0.0_r8

      SELECT CASE (fileformat)
       CASE DEFAULT
         avg_kernels_1(1:tropomi_dim)  = read_tropomi_avg_kernels(ifile, tropomi_dim, fileformat)
         obs_p_1(1:tropomi_dim) = read_tropomi_pressure(ifile, tropomi_dim, fileformat)
      END SELECT
      ! Dummy implementation, replace with actual code to read data from the file

      counts1 = counts1 + 1
      key = counts1
      call set_obs_def_tropomi(key, avg_kernels_1,obs_p_1)

   end subroutine read_tropomi_no2


   subroutine write_tropomi_no2(key, ifile, fform)
      integer, intent(in) :: key
      integer, intent(in)  :: ifile
      character(len=*),  intent(in), optional :: fform

      if ( .not. module_initialized ) call initialize_module
      if (ascii_file_format(fform)) then
         call write_tropomi_avg_kernels(ifile, kernel_trop_px, key, tropomi_dim, fform)
         write(ifile,11) key
         write(ifile, *) pressure_px(:,key)
11       format('pressure_px', i8)

      else
         call write_tropomi_avg_kernels(ifile, kernel_trop_px, key, tropomi_dim, fform)
         write(ifile,11) key
         write(ifile, *) pressure_px(:,key)
      end if
   end subroutine write_tropomi_no2

   subroutine set_obs_def_tropomi(key, kernel_trop_px_1, pressure_px_1)
      integer, intent(in) :: key
      integer :: i
      character(len=32)    :: fform
      real(r8),dimension(34) :: kernel_trop_px_1, pressure_px_1

      do i = 1, 34
         kernel_trop_px(i,key) = kernel_trop_px_1(i)
         pressure_px(i, key) = pressure_px_1(i)
      end do

   end subroutine set_obs_def_tropomi


!------------------------------
   subroutine ComputeWeight(nx,ny,xx,yy,ii,jj,dx,dy,nw,ww)
      implicit none
      real(r8),intent(in)        ::  xx(1:nx+1)  !(1:nx+1)
      real(r8),intent(in)        ::  yy(1:ny+1)  ! (1:ny+1)
      integer, intent(in)     ::  nx
      integer, intent(in)     ::  ny
      integer, intent(out)    :: ii(nx*ny)
      integer, intent(out)    :: jj(nx*ny)
      real(r8), intent(out)       :: ww(nx*ny)
      real(r8), intent(out)       :: dx(nx)
      real(r8), intent(out)       :: dy(ny)
      integer, intent(out)       :: nw


      ! --- local ----------------------------------

      integer     ::  xdir, ydir
      integer     ::  i, j
      integer     ::  is, in
      integer     ::  js, jn
      integer     ::  i1, i2
      real        ::  y1, y2
      real        ::  w
      real        ::  wsum

      ! --- begin ----------------------------------
      ii=0
      jj=0
      ww=0
      dx=0
      dy=0

      ! source cell range and step in increasing direction:
      if ( xx(1) < xx(nx+1) ) then
         is   = 1
         in   = nx
         xdir = 1
      else
         is   = nx
         in   = 1
         xdir = -1
      end if

      ! target cell range and step in increasing direction:
      if ( yy(1) < yy(ny+1) ) then
         js   = 1
         jn   = ny
         ydir = 1
      else
         js   = ny
         jn   = 1
         ydir = -1
      end if
      ! interval lengths:
      dx = xdir * ( xx(2:nx+1) - xx(1:nx) )
      ! check ..
      if ( any(dx <= 0.0) ) then
         write(*,'("found zero or negative interval length in x:")')
         do i = 1, nx
            write(*,'("  cell ",i6," [",f16.6,",",f16.6,"] ",f16.6)')       &
               i,xx(i), xx(i+1), dx(i)
         end do
         return
      end if


      ! interval lengths:
      dy = ydir * ( yy(2:ny+1) - yy(1:ny) )
      ! check ..
      if ( any(dy <= 0.0) ) then
         write(*,'("found zero or negative interval length in y:")')
         do i = 1, ny
            write(*,'("  cell ",i6," [",f16.6,",",f16.6,"] ",f16.6)')       &
               i,yy(i), yy(i+1), dy(i)
         end do
         return
      end if

      ! reset weights:
      nw = 0
      ! no first source cell yet:
      i1 = -999
      ! loop over target elements:
      do j = js, jn, ydir
         ! target cell bounds:
         if ( ydir > 0 ) then
            y1 = yy(j)
            y2 = yy(j+1)
         else
            y1 = yy(j+1)
            y2 = yy(j)
         end if

         ! search x cell holding y1 if not copied from i2:
         if ( i1 < 0 ) then
            ! loop over source cells until found:
            do i = is, in, xdir
               ! in cell bounds?
               if ( (xdir*xx(i) <= xdir*y1) .and.                            &
                  (xdir*y1 <= xdir*xx(i+1)) )  then
                  i1 = i
                  exit
               end if
            end do ! i
         end if ! search i1
         ! check ..
         if ( i1 < 0 ) then
            write(*,'("could not find interval holding y1 = ",f16.6)')y1
            do i = 1, nx
               write(*,'("  cell ",i6," [",f16.6,",",f16.6,"]")')            &
                  i, xx(i), xx(i+1)
            end do
            return
         end if

         ! search x cell holding y2:
         i2 = -999
         do i = i1, in, xdir
            if ( (xdir*xx(i) <= xdir*y2) .and.                              &
               (xdir*y2 <= xdir*xx(i+1)) ) then
               i2 = i
               exit
            end if
         end do ! i
         ! check ..
         if ( i2 < 0 ) then
            write(*,'("could not find interval holding y2 = ",f16.6)') y2
            do i = 1, nx
               write(*,'("  cell ",i6," [",f16.6,",",f16.6,"]")')             &
                  i,xx(i), xx(i+1)
            end do
            return
         end if

         do i = i1, i2, xdir
            ! fraction of source cell that overlaps with [y1,y2]:
            if ( xdir > 0 ) then
               w = ( min(xx(i+1),y2) - max(xx(i  ),y1) )/dx(i)
            else
               w = ( min(xx(i  ),y2) - max(xx(i+1),y1) )/dx(i)
            end if
            !! testing ..
            ! increase counter:
            nw = nw + 1
            ! store:
            ii(nw) = i
            jj(nw) = j
            ww(nw) = w
         end do ! i

         ! next:
         i1 = i2
      end do ! j

      ! check ..
      wsum = sum(ww(1:nw))
      if ( abs(wsum - 1.0) <= 1.0e-4 ) then
         write(*,'("sum of weights is ",es16.6," instead of 1.0")')wsum
         return
      end if
   end subroutine ComputeWeight

   subroutine FindIndexSurfacePressure(nz,pressure,indexSP)
      implicit none
      integer,intent(in) :: nz
      real(r8),intent(in) :: pressure(nz)
      integer,intent(out) :: indexSP
      if ( pressure(1) > pressure(nz) ) then
         indexSP = 1
      else
         indexSP = nz
      end if
   end subroutine FindIndexSurfacePressure

   function read_tropomi_avg_kernels(ifile, nlevels, fform)

      integer,                    intent(in) :: ifile, nlevels
      real(r8), dimension(tropomi_dim)        :: read_tropomi_avg_kernels
      character(len=*), intent(in), optional :: fform
      integer :: keyin
      character(len=14)   :: header
      character(len=129) :: errstring
      character(len=32)  :: fileformat

      read_tropomi_avg_kernels(:) = 0.0_r8

      if ( .not. module_initialized ) call initialize_module

      fileformat = "ascii"    ! supply default
      if(present(fform)) fileformat = trim(adjustl(fform))


      read(ifile, FMT='(a14, i8)') header, keyin    ! throw away keyin
      if(header /= 'kernel_trop_px') then
         call error_handler(E_ERR,'read tropomi avg kernels', &
            'Expected header "kernel_trop_px" in input file', source, revision, revdate)
      endif
      read(ifile, *) read_tropomi_avg_kernels(1:nlevels)

   end function read_tropomi_avg_kernels


   function read_tropomi_pressure(ifile, nlevels, fform)

      integer,                    intent(in) :: ifile, nlevels
      real(r8), dimension(tropomi_dim)        :: read_tropomi_pressure
      character(len=*), intent(in), optional :: fform

      character(len=11)   :: header
      character(len=129) :: errstring
      character(len=32)  :: fileformat
      integer :: keyin
      read_tropomi_pressure(:) = 0.0_r8

      if ( .not. module_initialized ) call initialize_module

      fileformat = "ascii"    ! supply default
      if(present(fform)) fileformat = trim(adjustl(fform))



      read(ifile, FMT='(a11, i8)') header, keyin    ! throw away keyin
      if(header /= 'pressure_px') then
         call error_handler(E_ERR,'read tropomi pressure_px', &
            'Expected header "pressure_px" in input file', source, revision, revdate)
         read(ifile, *) read_tropomi_pressure(1:nlevels)
      else
         read(ifile, *) read_tropomi_pressure(1:nlevels)
      end if        ! read and throw away

   end function read_tropomi_pressure



   subroutine write_tropomi_avg_kernels(ifile, avg_kernels_temp, key, nlevels_temp, fform)

      integer,                    intent(in) :: ifile, nlevels_temp
      real(r8), dimension(tropomi_dim), intent(in)  :: avg_kernels_temp
      character(len=32),          intent(in) :: fform

      character(len=5)   :: header
      character(len=129) :: errstring
      character(len=32)  :: fileformat
      integer :: key

      if ( .not. module_initialized ) call initialize_module

      fileformat = trim(adjustl(fform))

      if ( .not. module_initialized ) call initialize_module
      if (ascii_file_format(fform)) then
         write(ifile,11) key
         write(ifile, *) kernel_trop_px(:,key)
11       format('kernel_trop_px', i8)


      else
         write(ifile,11) key
         write(ifile, *) kernel_trop_px(:,key)
      end if

   end subroutine write_tropomi_avg_kernels


end module obs_def_SAT_NO2_TROPOMI_mod

! END DART PREPROCESS MODULE CODE

