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
! kernel and pressure vectors via set_obs_def_no2_tropomi. It subsequently calls 'write_tropomi_no2'
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
!   set_obs_def_no2_tropomi
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

! BEGIN DART PREPROCESS set_obs_def_no2_tropomi
!      case(SAT_NO2_TROPOMI)
!         call set_obs_def_no2_tropomi(obs_def%key)
! END DART PREPROCESS set_obs_def_no2_tropomi

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

   public ::  get_expected_SAT_NO2_TROPOMI, write_tropomi_no2,read_tropomi_no2, set_obs_def_no2_tropomi

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
   integer, parameter :: nretr = 1 ! hardcoded
   integer, parameter, dimension(nretr) :: nla  = 34! hardcoded
   integer, parameter :: max_obs = 1000000 ! number of intervals if model_levels is F
   integer,  dimension(max_obs)   :: tropomi_nlevels
   real(r8),dimension(tropomi_dim, max_obs) :: kernel_trop_px
   real(r8),dimension(tropomi_dim, max_obs) :: pressure_px
   real(r8),dimension(max_obs) :: amf
   character(len=6), parameter :: S5Pstring = 'FO_params'
   logical :: amf_from_model = .true.
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
      ! Input parameters
      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: ens_size
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: key
      ! Output parameters
      real(r8),            intent(out) :: val(ens_size)
      integer,             intent(out) :: istatus(ens_size)

      ! Local variables
      integer                          :: imem, num_levs, nz, iz, level_ith
      real(r8)                         :: model_p(ens_size, max_model_p_levs)
      real(r8)                         :: model_conc(ens_size, max_model_levs)
      real(r8)                         :: model_conc_2d_kl(ens_size, tropomi_dim)
      real(r8)                         :: Sx(ens_size, nretr)
      real(r8)                         :: model_conc_vcd(ens_size)
      real(r8)                         :: model_conc_vcd_alt(ens_size)
      real(r8)                         :: tropomi_pres_local(ens_size,tropomi_dim +1)
      real(r8)                         :: tropomi_trop_kernel_local(ens_size,tropomi_dim)
      real(r8)                         :: amf_local(ens_size)
      integer                          :: p_col_istatus(ens_size), int_conc_status(ens_size)
      type(location_type)              :: locS
      real(r8)                         :: mloc(3)
      real(r8)                         :: sp(ens_size)
      integer                          :: indexSP_kernel,indexSP_model
      integer                          :: iiv(tropomi_dim*max_model_levs),jjv(tropomi_dim*max_model_levs),nwv
      real(r8)                         :: dxv(max_model_levs),dyv(max_model_levs),wwv(tropomi_dim*max_model_levs)
      ! !Vertical for LayerAverage
      ! logical                          :: return_now
      ! if ( .not. module_initialized ) call initialize_module

      val = 0.0_r8
!     Get location information
      mloc = get_location(location)
      if (mloc(2)>90.0_r8) then
         mloc(2)=90.0_r8
      elseif (mloc(2)<-90.0_r8) then
         mloc(2)=-90.0_r8
      endif

      ! Initialize arrays
      tropomi_pres_local = 0.0_r8
      tropomi_trop_kernel_local = 0.0_r8

      ! Populate tropomi_pres_local and tropomi_trop_kernel_local arrays
      do imem = 1, ens_size
         do level_ith = 1, tropomi_dim
            tropomi_pres_local(imem, level_ith) = pressure_px(level_ith,key)
            tropomi_trop_kernel_local(imem, level_ith) = kernel_trop_px(level_ith, key)
         enddo
      enddo

      amf_local = amf(key)

      ! Initialize istatus
      istatus = 0
      nz = 1
      !     FARM pressure field at pixel position
      model_p = 0.0_r8
      model_levels: do
         if(nz == max_model_p_levs) then
            istatus = 0
            exit model_levels
         end if
         locS = set_location(mloc(1),mloc(2),farm_heights(nz),VERTISHEIGHT)
         call interpolate(state_handle, ens_size, locS, QTY_PRESSURE, model_p(:, nz), p_col_istatus)
         if (any(p_col_istatus /= 0)) return
         nz = nz + 1
      enddo model_levels
      ! nz should be at this point 16
      ! num_levs u
!     FARM surface pressure field at pixel position
      sp = 0.0_r8
      call interpolate(state_handle, ens_size, locS, QTY_SURFACE_PRESSURE, sp(:), p_col_istatus)
      if (any(p_col_istatus /= 0)) then
         sp(:) = MISSING_R8
      endif

      model_conc = 0.0_r8
      nz = 1
      !     FARM species field at pixel position
      model_levels_conc: do
         locS = set_location(mloc(1),mloc(2),farm_heights(nz),VERTISHEIGHT)
         call interpolate(state_handle, ens_size, locS, QTY_NO2, model_conc(:, nz), int_conc_status)
         ! call track_status(ens_size, int_conc_status, model_conc(:, nz), istatus, return_now)
         ! if(return_now) return
         if (any(int_conc_status /= 0)) then
            model_conc(:, nz) = MISSING_R8
            num_levs = nz - 1
            if(nz == 17) then
               istatus = 0
            end if
            exit model_levels_conc
         endif
         nz = nz + 1
      enddo model_levels_conc

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


      call FindIndexSurfacePressure(tropomi_dim,tropomi_pres_local(1, :), indexSP_kernel)
      call FindIndexSurfacePressure(nz+1,model_p(1, :),indexSP_model)

      do imem = 1, ens_size
         model_p(imem, :) = model_p(imem, :) / model_p(imem,indexSP_model) * tropomi_pres_local(imem, indexSP_kernel)
         model_p(imem, 1)=model_p(imem, 1)+0.1
      enddo
      ! Compute weights for vertical interpolation
      iiv = 0
      jjv = 0
      dxv = 0.0_r8
      dyv = 0.0_r8
      nwv = 0
      wwv  = 0.0_r8
      call ComputeWeight(max_model_levs,tropomi_dim,model_p(1, :),tropomi_pres_local(1,:),                   & !input
         iiv(:),jjv(:),dxv(:),dyv(:),nwv,wwv(:))
      model_conc_2d_kl = 0.0_r8
      do imem = 1, ens_size
         call ApplyWeightedSum(max_model_levs, tropomi_dim, nwv, iiv(:), jjv(:), wwv(:), dxv(:), &
            model_conc(imem,:), model_conc_2d_kl(imem, :))
      enddo

      ! conversion ug/m3 * Pa to mol/m2
      do imem = 1, ens_size
         do level_ith = 1, tropomi_dim
            model_conc_2d_kl(imem, level_ith) = UnitConversion_ugm3Pa_molm2(model_conc_2d_kl(imem, level_ith))
         end do
      end do
      model_conc_vcd = 0.0_r8
      model_conc_vcd_alt = 0.0_r8
      do imem = 1, ens_size
         call ApplyKernel(tropomi_dim, tropomi_trop_kernel_local(imem, :), model_conc_2d_kl(imem, :), model_conc_vcd(imem))
      end do
      ! vcd using the standard air mass factor
      do imem = 1, ens_size
         ! compute partial column sx
         call PartialColumn(nretr,tropomi_dim,nla,model_conc_2d_kl(imem, :),Sx(imem,:))
         ! compute alternative air mass factor correction M_m
         call AltAirMassFactor(nretr,tropomi_dim,amf_local,tropomi_trop_kernel_local(imem, :),model_conc_2d_kl(imem, :),Sx(imem,:),model_conc_vcd_alt(imem))
      end do

      val = model_conc_vcd
      istatus = 0

   end subroutine get_expected_SAT_NO2_TROPOMI

   subroutine read_tropomi_no2(key, ifile, fform)
      integer, intent(out) :: key
      integer :: i
      integer, intent(in)  :: ifile
      character(len=*), intent(in), optional    :: fform
      character(len=32) 		:: fileformat
      real(r8), dimension(tropomi_dim)	:: avg_kernels_1
      real(r8), dimension(tropomi_dim)	:: obs_p_1
      real(r8) :: amf_1
      integer 			:: keyin

      if ( .not. module_initialized ) call initialize_module

      fileformat = "ascii"   ! supply default
      if(present(fform)) fileformat = trim(adjustl(fform))


      avg_kernels_1(:) = 0.0_r8

      SELECT CASE (fileformat)
       CASE DEFAULT
         avg_kernels_1(1:tropomi_dim)  = read_tropomi_avg_kernels(ifile, tropomi_dim, fileformat)
         obs_p_1(1:tropomi_dim) = read_tropomi_pressure(ifile, tropomi_dim, fileformat)
         amf_1 = read_tropomi_amf(ifile, fform=fform)
      END SELECT
      ! Dummy implementation, replace with actual code to read data from the file

      counts1 = counts1 + 1
      key = counts1
      call set_obs_def_no2_tropomi(key, avg_kernels_1,obs_p_1, amf_1)

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
         call write_tropomi_amf(ifile, key, fform)
      else
         call write_tropomi_avg_kernels(ifile, kernel_trop_px, key, tropomi_dim, fform)
         write(ifile,11) key
         write(ifile, *) pressure_px(:,key)
         write(ifile, *) amf(key)
      end if
   end subroutine write_tropomi_no2

   subroutine set_obs_def_no2_tropomi(key, kernel_trop_px_1, pressure_px_1, amf_1)
      integer, intent(in) :: key
      integer :: i
      character(len=32)    :: fform
      real(r8),dimension(tropomi_dim) :: kernel_trop_px_1, pressure_px_1
      real(r8) :: amf_1

      do i = 1, tropomi_dim
         kernel_trop_px(i,key) = kernel_trop_px_1(i)
         pressure_px(i, key) = pressure_px_1(i)
         amf(key) = amf_1
      end do

   end subroutine set_obs_def_no2_tropomi


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

   function read_tropomi_amf(ifile, fform)

      integer,                    intent(in) :: ifile
      real(r8)        :: read_tropomi_amf
      character(len=*), intent(in), optional :: fform
      integer :: keyin
      character(len=14)   :: header
      character(len=129) :: errstring
      character(len=32)  :: fileformat

      read_tropomi_amf = 0.0_r8

      if ( .not. module_initialized ) call initialize_module

      fileformat = "ascii"    ! supply default
      if(present(fform)) fileformat = trim(adjustl(fform))


      read(ifile, FMT='(a7, i8)') header, keyin    ! throw away keyin
      if(header /= 'amf_tm5') then
         call error_handler(E_ERR,'read tropomi amf', &
            'Expected header "amf_tm5" in input file', source, revision, revdate)
      endif
      read(ifile, *) read_tropomi_amf

   end function read_tropomi_amf

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

   subroutine write_tropomi_amf(ifile, key, fform)

      integer,                    intent(in) :: ifile
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
         write(ifile, *) amf(key)
11       format('amf_tm5', i8)


      else
         write(ifile,11) key
         write(ifile, *) amf(key)
      end if

   end subroutine write_tropomi_amf


!-----------------------------------------------------------------
   subroutine ApplyWeightedSum(nx,ny,nw,ii, jj, ww, dx, f, g )
      !
      ! Compute partial sums:
      !
      !   g(j) = sum f(i) * w(i,j)
      !           i
      !
      ! with w(i,j) the faction of source interval i
      ! that overlaps with target interval j.
      !
      implicit none

      ! --- in/out ---------------------------------

      integer, intent(in)       ::  nx
      integer, intent(in)       ::  ny
      integer, intent(in)       ::  nw
      integer, intent(in)       ::  jj(nx*ny)
      integer, intent(in)       ::  ii(nx*ny)
      real(r8), intent(in)       ::  dx(nx)
      real(r8), intent(in)       ::  ww(nx*ny)
      real(r8), intent(in)       ::  f(nx)  ! (1:nx)
      real(r8), intent(out)      ::  g(ny)  ! (1:ny)

      ! --- const ----------------------------------


      ! --- local ----------------------------------

      integer     ::  iw

      ! --- begin ----------------------------------


      ! init result:
      g = 0.0
      ! loop over mapping weights:
      do iw = 1, nw
         ! add contribution:
         g(jj(iw)) = g(jj(iw)) + f(ii(iw)) * dx(ii(iw))                  &
            * ww(iw)
!        g(jj(iw)) = g(jj(iw)) + f(ii(iw)) * ww(iw)

      end do ! iw

   end subroutine ApplyWeightedSum


!------------------------------------------------------------
   real function UnitConversion_ppbPa_molm2(varin)

      !if varin = ppb*Pa -> output of ApplyWeightedSum to concentration in ppb then
      !  (mole tr)/m2 =
      ! ppb*Pa
      ! / [g]   :  Pa/[g] = (kg air/m2)
      ! / ((kg air)/(mole air))
      ! * (mole tr)/(mole  air)/ppb

      !if varin =  -> output of ApplyWeightedSumAdj to gradient (xo-Hyb)/(HBHT+R) then
      ! 1/ppb =
      ! Pa/[(mol tr)/m2]
      ! / [g]   :  Pa/[g] = (kg air)/m2
      ! / ((kg air)/(mole air))
      ! * (mole tr)/(mole * air)/ppb
      implicit none
      ! in/out --------------------------------
      real     :: varin
      ! local ----------------------------------
      ! gravity constant:
      real, parameter     ::  grav   = 9.80665 ! m/s2
      ! mole mass of air:
      real, parameter     ::  xm_air = 28.964e-3     ! kg/mol : ~80% N2, ~20% O2
      ! begin ----------------------------------

      ! unit conversion:
      UnitConversion_ppbPa_molm2         =   &   !  (mole tr)/m2 =
         varin   &   ! ppb*Pa
         / grav   &   ! / [g]   :  Pa/[g] = (kg air/m2)
         / xm_air &   ! / ((kg air)/(mole air))
         * 1.0e-9      ! * (mole tr)/(mole  air)/ppb
   end function UnitConversion_ppbPa_molm2


   real function UnitConversion_ugm3Pa_molm2(varin)

      implicit none
      ! in/out --------------------------------
      real(r8)     :: varin
      ! local ----------------------------------
      ! gravity constant:
      real, parameter     ::  grav   = 9.80665 ! m/s2
      ! mole mass of air:
      real, parameter     ::  rho_air = 1.225     ! kg_air/m3
      real, parameter     ::  Mw = 46.0055     ! g_tr/mol_tr
      ! begin ----------------------------------

      ! unit conversion:
      UnitConversion_ugm3Pa_molm2         =   &   !  (mole tr)/m2 =
         varin   &   ! ug/m3*kg/m*s2
         / grav   &   ! m/s2
         / rho_air & ! kg/m3
         / Mw &   ! / ((g air)/(mole air))
         * 1.0e-6      ! * g / ug
   end function UnitConversion_ugm3Pa_molm2


!---------------------------------------------------------
   subroutine ApplyKernel(nlayer,A_data,x_data,y_data)
      implicit none
      integer,intent(in)   :: nlayer
      integer :: layeri
      real(r8), intent(in)     :: A_data(nlayer)
      real(r8), intent(in)     :: x_data(nlayer)
      real(r8), intent(out)    :: y_data

      ! --- begin ----------------------------------
      y_data = 0.0_r8
      do layeri = 1, nlayer
         y_data = y_data + A_data(layeri)*x_data(layeri)
      end do
   end subroutine ApplyKernel
!---------------------------------------------------------

   subroutine PartialColumn(nretr,nlayer,nla,x,Sx)

      !-in/out ----------------
      integer,intent(in)     :: nlayer
      integer,intent(in)     :: nretr
      integer,intent(in)     :: nla(nretr)
      real(r8),intent(in)        :: x(nlayer)
      real(r8),intent(out)       :: Sx(nretr)
      ! local -----------------
      integer                :: k1,k2,iretr
      ! begin -----------------
      ! init index:
      k2 = 0
      ! loop over retrieval layers
      do iretr = 1, nretr
         ! range of apriori layers
         k1 = k2 + 1
         k2 = k2 + nla(iretr)
         ! fill:
         Sx(iretr) = sum(x(k1:k2))
      end do ! iretr
   end subroutine PartialColumn
   !---------------------------------------------------------
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !~ airmass factor from local model:
   !    M_m = M A x / (sum_{l=1,nla} x)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine AltAirMassFactor(nretr,nlayer,M, A, x, Sx, M_m)
      implicit none
      ! in/out -----------------------------------
      integer,intent(in)    :: nretr
      integer,intent(in)    :: nlayer
      real(r8),intent(in)       :: M(nretr)
      real(r8),intent(in)       :: A(nretr,nlayer)
      real(r8),intent(in)       :: x(nlayer)
      real(r8),intent(in)       :: Sx(nretr)
      real(r8),intent(out)      :: M_m(nretr)
      ! local -----------------------------------
      integer    :: iretr
      ! begin -----------------------------------

      ! loop over retrieval layers
      do iretr = 1, nretr
         ! fill, trap divison by zero:
         if ( Sx(iretr) >= 0.0 ) then
            M_m(iretr) = M(iretr) *      &
               sum(A(iretr,:) * x(:) )    &
               / Sx(iretr)
         else
            M_m(iretr) = M(iretr)
         end if
      end do ! iretr
   end subroutine AltAirMassFactor
   !------------------------------------------------------------
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !~ averaging kernel using airmass factor from local model:
   !    A_m = M / M_m  A
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine AltKernel(nretr, nlayer, A, M, M_m, A_m)
      implicit none
      ! in/out ------------------------
      integer, intent(in) :: nretr
      integer, intent(in) :: nlayer
      real, intent(in)    :: A(nretr,nlayer)
      real, intent(in)    :: M(nretr)
      real, intent(in)    :: M_m(nretr)
      real, intent(out)    :: A_m(nretr,nlayer)
      ! local -----------------------------------
      integer    :: iretr
      ! begin -----------------------------------
      ! loop over retrieval layers
      do iretr = 1, nretr
         if (  M_m(iretr) >= 0.0 )   then
            A_m(iretr,:) = M(iretr) / M_m(iretr) * A(iretr,:)
         else
            A_m(iretr,:) =  A(iretr,:)
         end if
      enddo
   end subroutine AltKernel
   !------------------------------------------------------------
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !~ retrieval using airmass factor from local model:
   !    y_m = M / M_m  y
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine AltRetrieval(nretr, y, M, M_m, y_m )
      implicit none
      ! in/out ---------------------------
      integer, intent(in) :: nretr
      real, intent(in) :: y(nretr)
      real, intent(in) :: M(nretr)
      real, intent(in) :: M_m(nretr)
      real, intent(out) :: y_m(nretr)
      ! local  ---------------------------
      integer    :: iretr
      ! begin ---------------------------
      ! loop over retrieval layers
      do iretr = 1, nretr
         if ( M_m(iretr) >= 0.0 ) then
            y_m(iretr) = M(iretr) / M_m(iretr) *y(iretr)
         else
            y_m(iretr) =  y(iretr)
         endif
      end do ! iretr
   end subroutine AltRetrieval
   !------------------------------------------------------------
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !~ retrieval error covariance using airmass factor from local
   !model:
   !    R_m = M / M_m  R
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine AltRetrievalCovar(nretr0,nretr, R, M, M_m, R_m )
      implicit none
      ! in/out -------------------------
      integer, intent(in) :: nretr0
      integer, intent(in) :: nretr
      real, intent(in)    :: R(nretr0,nretr)
      real, intent(in)    :: M(nretr)
      real, intent(in)    :: M_m(nretr)
      real, intent(out)   :: R_m(nretr0,nretr)
      ! local begin --------------------
      integer :: k1,k2
      ! begin --------------------------
      do k2 = 1, nretr
         do k1 = 1, nretr0
            if ( ( M_m(k1) >= 0.0 ) .and. M_m(k2) >= 0 ) then
               R_m(k1,k2) = M(k1)/M_m(k1) * R(k1,k2) * M(k2)/M_m(k2)
            else
               R_m(k1,k2) =   R(k1,k2)
            endif
         end do !iretr0
      end do ! iretr
   end subroutine AltRetrievalCovar
   !------------------------------------------------------------
   real function UnitConversion(varin)

      !if varin = ppb*Pa -> output of ApplyWeightedSum to concentration in ppb then
      !  (mole tr)/m2 =
      ! ppb*Pa
      ! / [g]   :  Pa/[g] = (kg air/m2)
      ! / ((kg air)/(mole air))
      ! * (mole tr)/(mole  air)/ppb

      !if varin =  -> output of ApplyWeightedSumAdj to gradient (xo-Hyb)/(HBHT+R) then
      ! 1/ppb =
      ! Pa/[(mol tr)/m2]
      ! / [g]   :  Pa/[g] = (kg air)/m2
      ! / ((kg air)/(mole air))
      ! * (mole tr)/(mole * air)/ppb
      implicit none
      ! in/out --------------------------------
      real     :: varin
      ! local ----------------------------------
      ! gravity constant:
      real, parameter     ::  grav   = 9.80665 ! m/s2
      ! mole mass of air:
      real, parameter     ::  xm_air = 28.964e-3     ! kg/mol : ~80% N2, ~20% O2
      ! begin ----------------------------------

      ! unit conversion:
      UnitConversion         =   &   !  (mole tr)/m2 =
         varin   &   ! ppb*Pa
         / grav   &   ! / [g]   :  Pa/[g] = (kg air/m2)
         / xm_air &   ! / ((kg air)/(mole air))
         * 1.0e-9      ! * (mole tr)/(mole  air)/ppb
   end function UnitConversion
   !------------------------------------------------------------
   real*8 function distlola(lon1,lat1,lon2,lat2)
      !horizontal distance (km) from lon,lat
      implicit none
      ! in/out ----------------
      real*8  :: lon1,lat1,lon2,lat2
      ! local  ----------------
      real*8,parameter :: RT = 6373.0447000
      real*8,parameter :: rtd = dasin(dble(1.)/dble(90.))
      real*8           :: rlat1,rlat2,rlon1,rlon2
      real*8           :: gam,alpha
      ! begin --------------------
      if ( ( lat1 < -90.0 ) .or. ( lat1 > 90.0 )  .or. &
         ( lat2 < -90.0 ) .or. ( lat2 > 90.0 )  .or. &
         ( lon1 < -360.0) .or. ( lon1 > 360.0)  .or. &
         ( lon2 < -360.0) .or. ( lon2 > 360.0)  ) then
         distlola = -999
      else if ( ( lat2 .eq. lat1) .and. ( lon2 .eq. lon1) ) then
         distlola = 0.0
      else
         rlat1 = rtd *lat1
         rlat2 = rtd *lat2
         rlon1 = rtd *lon1
         rlon2 = rtd *lon2
         gam = dcos(rlat1) *dcos(rlat2) *dcos(rlon1 -rlon2) + dsin(rlat1) *dsin(rlat2)
         alpha = dacos(gam)
         distlola = RT* alpha
      endif
   end function distlola
   !------------------------------------------------------------
   subroutine BilinearInt(nz,nxc,nyc,nxf,nyf,           &!dimensions
      fieldc,idlon,idlat,wlon,wlat, &!input
      fieldf)                        !output
      implicit none
      ! -- i/o var ---
      integer :: nxc,nyc,nxf,nyf,nz
      real    :: fieldc(nxc,nyc,nz)
      real    :: fieldf(nxf,nyf,nz)
      integer :: idlon(2,nxf),idlat(2,nyf)
      real    :: wlon(nxf),wlat(nyf)
      ! -- local ----
      integer :: i,j,k

      ! -- begin ----
      fieldf(:,:,:)=0.0
      !$OMP PARALLEL DO                                     &
      !$OMP COLLAPSE(3)                                     &
      !$OMP DEFAULT (NONE)                                  &
      !$OMP PRIVATE(i,j,k)                                  &
      !$OMP SHARED(nz,nyf,nxf)                              &
      !$OMP SHARED(fieldf,fieldc,idlon,idlat,wlon,wlat)
      do k=1,nz
         do j=1,nyf
            do i=1,nxf
               fieldf(i,j,k)                                          = &
                  fieldc(idlon(1,i),idlat(1,j),k)*wlon(i)*wlat(j)        + &
                  fieldc(idlon(2,i),idlat(2,j),k)*(1-wlon(i))*(1-wlat(j))+ &
                  fieldc(idlon(1,i),idlat(2,j),k)*(wlon(i))*(1-wlat(j))  + &
                  fieldc(idlon(2,i),idlat(1,j),k)*(1-wlon(i))*(wlat(j))
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP BARRIER

   end subroutine BilinearInt
   !------------------------------------------------------------
   subroutine sorting_dble(n,array,sort_index)
      implicit none
      ! --- i/o vars ----
      integer,intent(in) :: n
      real*8,intent(in)  :: array(n)
      integer,intent(out) :: sort_index(n)
      ! --- local ---
      integer :: i,j,m,idx
      real*8  :: mindist
      real*8  :: array2(n)
      ! --- begin ----
      array2 = array
      mindist =  -999
      idx=1
      do i=1,n
         mindist = minval(array2,mask=(array2>mindist))
         m = count(array2==mindist)
         do j = 1,m
            sort_index(idx) = minloc(array2,dim=1,mask=array2==mindist)
            array2(sort_index(idx))=-999
            idx = idx+1
         enddo
      enddo
   end subroutine sorting_dble

end module obs_def_SAT_NO2_TROPOMI_mod

! END DART PREPROCESS MODULE CODE

