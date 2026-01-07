module nc_interface_mod

contains
   !---------------------------------------------
   subroutine Nc_Getattr(ncid,varname,attrname,multiplication_factor)
      use NetCDF , only :  NF90_inq_varid, NF90_Inquire_Variable,       &
         NF90_inquire_dimension,NF90_get_var,         &
         NF90_get_att
      IMPLICIT NONE
      real(kind=kind(0.0)) :: multiplication_factor
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: attrname
      integer                     :: ierr
      integer,intent(in)          :: ncid
      integer                     :: varid
      integer                     :: ndims,status
      integer,allocatable,dimension(:) :: dimlen
      integer,allocatable,dimension(:) :: dimids
      integer                     :: nc_status
      integer                             :: i

      ! Retrieve variable ID
      status = NF90_inq_varid(ncid, varname, varid)
      call handle_err(status)

      ierr = NF90_get_att(ncid, varid, 'multiplication_factor_to_convert_to_molecules_percm2', multiplication_factor)

   end subroutine Nc_Getattr
   !---------------------------------------------
   subroutine Nc_GetDim(ncid,dimname,length)

      use NetCDF , only :  NF90_Inq_DimID, NF90_Inquire_Dimension, NF90_get_att

      IMPLICIT NONE

      integer,intent(in)          :: ncid
      character(len=*),intent(in) :: dimname
      integer,intent(out)         :: length
      integer                     :: dimid
      integer                     :: status

      status = NF90_Inq_DimID( ncid,trim(dimname), dimid )
      call handle_err(status)
      status = NF90_Inquire_Dimension( ncid, dimid, len=length )
      call handle_err(status)

   end subroutine Nc_GetDim
   !---------------------------------------------

   subroutine Nc_GetVar(ncid,varname,var1d,var2d,var3d,var4d,ivar1d,ivar2d,timestep,neofs)

      use NetCDF , only :  NF90_inq_varid, NF90_Inquire_Variable,       &
         NF90_inquire_dimension,NF90_get_var,         &
         NF90_get_att
      IMPLICIT NONE
      real(kind=kind(0.0)) :: scale_factor, add_offset
      character(len=*),intent(in) :: varname
      integer                     :: ierr
      integer,intent(in)          :: ncid
      integer                     :: varid
      integer                     :: ndims,natts
      integer,allocatable,dimension(:) :: dimlen
      integer,allocatable,dimension(:) :: dimids
      integer                     :: status
      integer,optional,intent(in)  :: timestep
      integer,optional,intent(in)  :: neofs
      real,allocatable,dimension(:),optional,intent(out) :: var1d
      integer*8,allocatable,dimension(:),optional,intent(out) :: ivar1d
      real,allocatable,dimension(:,:),optional,intent(out) :: var2d
      integer,allocatable,dimension(:,:),optional,intent(out) :: ivar2d
      real,allocatable,dimension(:,:,:),optional,intent(out) :: var3d
      real,allocatable,dimension(:,:,:,:),optional,intent(out) :: var4d
      integer                             :: i

      status = NF90_inq_varid(ncid, varname, varid)
      call handle_err(status)

      ierr = NF90_get_att(ncid, varid, 'scale_factor', scale_factor)
      ierr = NF90_get_att(ncid, varid, 'add_offset', add_offset)

      if (scale_factor == 0.0) then
         scale_factor = 1.0
      end if

      status = NF90_Inquire_Variable( ncid, varid,                      &
         ndims=ndims,natts=natts)
      call handle_err(status)

      allocate(dimlen(ndims),dimids(ndims))

      status = NF90_Inquire_Variable( ncid, varid,dimids=dimids)
      call handle_err(status)

      do i = 1,ndims
         status = NF90_inquire_dimension(ncid,dimids(i),len = dimlen(i))
         call handle_err(status)
      enddo

      if ( ndims == 1 ) then
         if (present(var1d)) then
            if (present(neofs)) then
               if ( .not. allocated(var1d)) allocate(var1d(neofs))
               status = NF90_get_var(ncid, varid, var1d,start=(/1/),count=(/neofs/))
               call handle_err(status)
               var1d = var1d * scale_factor + add_offset
            else
               if ( .not. allocated(var1d)) allocate(var1d(dimlen(1)))
               status = NF90_get_var(ncid, varid, var1d)
               call handle_err(status)
               var1d = var1d * scale_factor + add_offset
            endif
         elseif (present(ivar1d) ) then
            if ( .not. allocated(ivar1d)) allocate(ivar1d(dimlen(1)))
            status = NF90_get_var(ncid, varid, ivar1d)
            call handle_err(status)
            ivar1d = ivar1d * scale_factor + add_offset
         endif
      elseif ( ndims == 2 ) then
         if (present(var2d)) then
            if (present(neofs)) then
               if ( .not. allocated(var2d)) allocate(var2d(dimlen(1),neofs))
               status = NF90_get_var(ncid, varid, var2d,start=(/1,1/),count=(/dimlen(1),neofs/))
               call handle_err(status)
               var2d = var2d * scale_factor + add_offset
            else
               if ( .not. allocated(var2d)) allocate(var2d(dimlen(1),dimlen(2)))
               status = NF90_get_var(ncid, varid, var2d)
               call handle_err(status)
               var2d = var2d * scale_factor + add_offset
            endif
         elseif (present(ivar2d)) then
            if ( .not. allocated(ivar2d)) allocate(ivar2d(dimlen(1),dimlen(2)))
            status = NF90_get_var(ncid, varid, ivar2d)
            call handle_err(status)
            ivar2d = ivar2d * scale_factor + add_offset
         endif
      elseif ( ndims == 3 ) then
         if (present(timestep)) then
            if ( .not. allocated(var2d)) allocate(var2d(dimlen(1),dimlen(2)))
            status = NF90_get_var(ncid, varid,var2d,                   &
               start=(/1,1,timestep/),                &
               count=(/dimlen(1),dimlen(2),1/))
            call handle_err(status)
            var2d = var2d * scale_factor + add_offset
         else
            if ( .not. allocated(var3d)) allocate(var3d(dimlen(1),dimlen(2),dimlen(3)))
            status = NF90_get_var(ncid, varid, var3d)
            var3d = var3d * scale_factor + add_offset
            call handle_err(status)
         end if
      elseif ( ndims == 4 ) then
         if (present(timestep)) then
            if ( .not. allocated(var3d)) allocate(var3d(dimlen(1),dimlen(2),dimlen(3)))
            status = NF90_get_var(ncid, varid,var3d,                       &
               start=(/1,1,1,timestep/),                &
               count=(/dimlen(1),dimlen(2),dimlen(3),1/))
            call handle_err(status)
            var3d = var3d * scale_factor + add_offset
         elseif (present(neofs)) then
            !            if ( .not. allocated(var4d)) allocate(var4d(dimlen(1),dimlen(2),dimlen(3),neofs))
            !            status = NF90_get_var(ncid, varid,var4d,                &
            !                                  start=(/1,1,1,1/),                &
            !                                  count=(/dimlen(1),dimlen(2),dimlen(3),neofs/))
            !            call handle_err(status)
            if ( .not. allocated(var4d)) allocate(var4d(neofs,dimlen(2),dimlen(3),dimlen(4)))
            status = NF90_get_var(ncid, varid,var4d,                &
               start=(/1,1,1,1/),                &
               count=(/neofs,dimlen(2),dimlen(3),dimlen(4)/))
         else
            if ( .not. allocated(var4d)) allocate(var4d(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
            status = NF90_get_var(ncid, varid, var4d)
            call handle_err(status)
            var4d = var4d * scale_factor + add_offset
         endif
      else
         print*,'ERROR: number of dimensions'
      endif

      deallocate(dimlen,dimids)
   end subroutine Nc_GetVar

   !-----------------------------------------------------
   subroutine Nc_GetUnits(ncid,varname,Units)

      use NetCDF , only :  NF90_inq_varid, NF90_Inquire_Variable,       &
         NF90_inq_attname, NF90_get_att

      IMPLICIT NONE
      character(len=*),intent(in) :: varname
      character(len=256)          :: name
      integer,intent(in)          :: ncid
      integer                     :: varid
      integer                     :: ndims,natts
      integer                     :: status
      integer                     :: i
      character(len=256),intent(out)   :: Units

      status = NF90_inq_varid(ncid, varname, varid)
      call handle_err(status)
      status = NF90_Inquire_Variable( ncid, varid,                      &
         ndims=ndims,natts=natts )
      do i = 1,natts
         status = NF90_inq_attname(ncid, varid, attnum=i, name=name)
         call handle_err(status)
         if ( name .ne. 'units' ) cycle
         status = NF90_get_att(ncid, varid, name, Units)
         call handle_err(status)
      enddo

   end subroutine Nc_GetUnits
   !-----------------------------------------------------

   subroutine Nc_open(filename,ncid)

      use NetCDF , only :  NF90_Open, NF90_NOWRITE

      IMPLICIT NONE

      character(len=*),intent(in) :: filename
      integer,intent(out)         :: ncid
      integer                     :: status

      status = NF90_Open(trim(filename), NF90_NOWRITE, ncid )
      call handle_err(status)

   end subroutine Nc_open
   !-----------------------------------------------------

   subroutine Nc_close(ncid)

      use NetCDF , only :  NF90_Close

      IMPLICIT NONE

      integer,intent(in)         :: ncid
      integer                     :: status

      status = NF90_close( ncid )
      call handle_err(status)

   end subroutine Nc_close
   !-----------------------------------------------------

   subroutine handle_err(status)

      use NetCDF , only : nf90_noerr, nf90_strerror

      IMPLICIT NONE

      integer, intent ( in) :: status
      if(status /= nf90_noerr) then
         print*,trim(nf90_strerror(status))
         stop "Stopped"
      end if

   end subroutine handle_err

end module nc_interface_mod

module DateTimeModule
   implicit none

   type :: DateTime
      integer :: year
      integer :: month
      integer :: day
      integer :: hour
      integer :: minute
      integer :: second
      integer :: millisecond
   end type DateTime

contains

   subroutine extractDateTime_from_timestr(timestamp_str, dt)
      character(len=*), intent(in) :: timestamp_str
      type(DateTime), intent(out) :: dt

      integer :: i
      character(len=4) :: year_str
      character(len=2) :: month_str, day_str, hour_str, minute_str, second_str, millisec_str

      ! Extract components from the timestamp string
      year_str = timestamp_str(20:23)
      month_str = timestamp_str(25:26)
      day_str = timestamp_str(28:29)
      hour_str = timestamp_str(31:32)
      minute_str = timestamp_str(34:35)
      second_str = timestamp_str(37:38)
      millisec_str = timestamp_str(40:48)

      ! Convert strings to integers
      read(year_str, *) dt%year
      read(month_str, *) dt%month
      read(day_str, *) dt%day
      read(hour_str, *) dt%hour
      read(minute_str, *) dt%minute
      read(second_str, *) dt%second
      read(millisec_str, *) dt%millisecond

   end subroutine extractDateTime_from_timestr

   subroutine extractDateTime_from_units(time, time_units, dt)
      implicit none

      ! Inputs
      integer(kind=8), intent(in) :: time
      character(len=*), intent(in) :: time_units

      ! Output
      type(DateTime), intent(out) :: dt

      ! Local
      integer :: ref_year, ref_month, ref_day
      integer :: ref_hour, ref_minute, ref_second
      integer :: days_in_month(12)
      integer(kind=8) :: remaining_seconds
      integer :: i
      integer :: ios

      character(len=32) :: dummy

      ! --------------------------------------------------
      ! Parse reference time from units string
      ! Expected: "seconds since YYYY-MM-DD HH:MM:SS"
      ! --------------------------------------------------

      ! Defaults if time part is missing
      ref_hour   = 0
      ref_minute = 0
      ref_second = 0
      dt%millisecond = 0

      ! NetCDF time units are assumed to follow the CF convention:
      !   "seconds since YYYY-MM-DD [HH:MM:SS[.mmm]]"
      !
      ! Date fields (Y-M-D) are mandatory.
      ! Time fields (H:M:S) are optional and default to zero if not present.
      ! Parsing is done via fixed-position slicing with iostat protection
      ! to avoid failures due to formatting variations.

      read(time_units(15:18), *, iostat=ios) ref_year
      if (ios /= 0) return

      read(time_units(20:21), *, iostat=ios) ref_month
      if (ios /= 0) return

      read(time_units(23:24), *, iostat=ios) ref_day
      if (ios /= 0) return

      read(time_units(48:49), *, iostat=ios) ref_hour
      if (ios /= 0) ref_hour = 0

      read(time_units(34:35), *, iostat=ios) ref_minute
      if (ios /= 0) ref_minute = 0

      read(time_units(37:38), *, iostat=ios) ref_second
      if (ios /= 0) ref_second = 0

      ! ---- INITIALIZE dt FROM REFERENCE ----
      dt%year        = ref_year
      dt%month       = ref_month
      dt%day         = ref_day
      dt%hour        = ref_hour
      dt%minute      = ref_minute
      dt%second      = ref_second

      remaining_seconds = time

      ! --------------------------------------------------
      ! Add seconds → hh:mm:ss
      ! --------------------------------------------------
      dt%second = dt%second + mod(remaining_seconds, 60)
      remaining_seconds = remaining_seconds / 60

      dt%minute = dt%minute + mod(remaining_seconds, 60)
      remaining_seconds = remaining_seconds / 60

      dt%hour = dt%hour + mod(remaining_seconds, 24)
      remaining_seconds = remaining_seconds / 24

      ! Normalize seconds/minutes/hours
      call normalize_hms(dt)

      ! --------------------------------------------------
      ! Add remaining days
      ! --------------------------------------------------
      do while (remaining_seconds > 0)

         if (is_leap_year(dt%year)) then
            days_in_month = (/31,29,31,30,31,30,31,31,30,31,30,31/)
         else
            days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         end if

         if (dt%day < days_in_month(dt%month)) then
            dt%day = dt%day + 1
         else
            dt%day = 1
            if (dt%month < 12) then
               dt%month = dt%month + 1
            else
               dt%month = 1
               dt%year = dt%year + 1
            end if
         end if

         remaining_seconds = remaining_seconds - 1
      end do

   end subroutine extractDateTime_from_units


   subroutine DateTimeInSeconds(dt)
      type(DateTime), intent(out) :: dt
      integer :: seconds
      seconds = 0
      ! Calculate total seconds
      seconds = dt%year * 365*24*3600 + dt%month * 30*24*3600 + dt%day * 24*3600 + &
         dt%hour * 3600 + dt%minute * 60 + dt%second
   end subroutine DateTimeInSeconds


   logical function is_leap_year(year)
      integer, intent(in) :: year
      is_leap_year = (mod(year,4) == 0 .and. mod(year,100) /= 0) .or. &
         (mod(year,400) == 0)
   end function is_leap_year

   subroutine normalize_hms(dt)
      type(DateTime), intent(inout) :: dt

      if (dt%second >= 60) then
         dt%minute = dt%minute + dt%second / 60
         dt%second = mod(dt%second, 60)
      end if

      if (dt%minute >= 60) then
         dt%hour = dt%hour + dt%minute / 60
         dt%minute = mod(dt%minute, 60)
      end if

      if (dt%hour >= 24) then
         dt%hour = mod(dt%hour, 24)
         ! days handled in main loop
      end if
   end subroutine normalize_hms


end module DateTimeModule
