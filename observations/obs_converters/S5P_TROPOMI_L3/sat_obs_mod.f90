module sat_obs_mod
   use DateTimeModule
   implicit none
   type :: T_SatObs
      ! number of pixels in global domain:
      integer :: npix
      ! number of corners in footprint:
      integer :: ncorner
      ! number of layers in profiles:
      integer ::  nlayer
      ! number of layers at interface:
      integer ::  nlayeri
      integer ::  nretr
      integer ::  nretr0
      !global coordinates of track
      real,allocatable :: lon(:)                   !(pixel)
      real,allocatable :: lat(:)                   !(pixel)
      real,allocatable :: clon(:,:)                !(ncorner,pixel)
      real,allocatable :: clat(:,:)                !(ncorner,pixel)
      real,allocatable :: qa_flag(:)               !(pixel)
      real,allocatable :: pressure(:,:)            !(layeri,pixel)
      real,allocatable :: kernel_trop(:,:,:)       !(retr,layer,pixel)
      real,allocatable :: amf_trop(:,:)            !(retr,pixel)
      real,allocatable :: vcd(:,:)                 !(retr,pixel)
      real,allocatable :: vcd_errvar(:,:,:)        !(retr0,retr,pixel)
      integer,allocatable :: nla(:,:)              !(retr,pixel)
      real,allocatable :: vcd_m(:,:)               !(retr,pixel)
      real,allocatable :: kernel_trop_m(:,:,:)     !(retr,layer,pixel)
      real,allocatable :: amf_trop_m(:,:)          !(retr,pixel)
      real,allocatable :: vcd_errvar_m(:,:,:)      !(retr0,retr,pixel)
      real,allocatable :: mdlspace_vcd_m(:,:,:)    !(nlon,nlat,1)
      real,allocatable :: mdlspace_vcd(:,:,:)      !(nlon,nlat,1)
      real,allocatable :: mdlspace_vr(:,:,:)       !(nlon,nlat,1)
      real,allocatable :: time(:)       !(nlon,nlat,1)
      type(DateTime), allocatable   :: date_time(:)
      integer,allocatable :: tobs(:) ! time in seconds
      character(len=256) :: time_units
   contains
      procedure ::  Read         =>  ReadSatObs
      procedure ::  Done         =>  SatObsDone
   end type

contains

   subroutine ReadSatObs(self,filename,pollutant)

      use nc_interface_mod , only : Nc_open,Nc_GetDim,  &
         Nc_GetVar,Nc_close, &
         Nc_GetUnits
      use DateTimeModule

      implicit none

      character(len=*),intent(in)       :: pollutant
      character(len=*),intent(in)       :: filename
      class(T_SatObs), intent(out)      :: self
      type(DateTime) :: reference_time, dt

!----------local
      integer                      :: ncid, i

!----------begin

!Open File
      call Nc_open(filename,ncid)
!Get Dimensions
      call Nc_GetDim(ncid,'pixel',self%npix)
      call Nc_GetDim(ncid,'corner',self%ncorner)
      call Nc_GetDim(ncid,'layer',self%nlayer)
      call Nc_GetDim(ncid,'layeri',self%nlayeri)
      call Nc_GetDim(ncid,'retr',self%nretr)
      call Nc_GetDim(ncid,'retr0',self%nretr0)
!Get Variables
      call Nc_GetVar(ncid,'longitude',var1d=self%lon)
      call Nc_GetVar(ncid,'latitude',var1d=self%lat)
      call Nc_GetVar(ncid,'time',var1d=self%time)
      call Nc_GetUnits(ncid, 'time', self%time_units)
      allocate(self%date_time(self%npix))
      do i = 1, size(self%time)
         ! A.D'A we should then add the time expressed in millisecond
         ! for the purposes of the assimilation in FARM the timestamp in the units is sufficient
         ! because data acquired within an hour belong to the same timestep.
         call extractDateTime(self%time_units, self%date_time(i))
         write(*, '(A, I4, A, I2, A, I2, A, I2, A, I2, A, I2, A, I2, A, I6)') &
            "Datetime: ", dt%year, "-", dt%month, "-", dt%day, " ", dt%hour, ":", dt%minute, ":", dt%second, ".", dt%millisecond
      end do
      call Nc_GetVar(ncid,'longitude_bounds',var2d=self%clon)
      call Nc_GetVar(ncid,'latitude_bounds',var2d=self%clat)
      call Nc_GetVar(ncid,'qa_value',var1d=self%qa_flag)
      call Nc_GetVar(ncid,'pressure',var2d=self%pressure)
      call Nc_GetVar(ncid,'vcd',var2d=self%vcd)
      call Nc_GetVar(ncid,'vcd_errvar',var3d=self%vcd_errvar)
      if ( pollutant == 'NO2' ) then
         call Nc_GetVar(ncid,'kernel_trop',var3d=self%kernel_trop)
         call Nc_GetVar(ncid,'amf_trop',var2d=self%amf_trop)
         call Nc_GetVar(ncid,'nla',ivar2d=self%nla)
      elseif ( ( pollutant == 'SO2' ) .or. ( pollutant == 'HCHO' ) ) then
         ! Kernel is defined up to the stratosphere
         ! Consistently with kernel and vcd amf is chosen for polluted scenario
         call Nc_GetVar(ncid,'kernel',var3d=self%kernel_trop)
         call Nc_GetVar(ncid,'amf',var2d=self%amf_trop)
         if ( .not. allocated(self%nla)) allocate(self%nla(self%nretr,self%npix))
         self%nla(:,:)=self%nlayer
      else
         print*,'WRONG pollutant: ',trim(pollutant)
      endif
!Close File
      call Nc_close(ncid)
!Teatment of Negative values
!     where ( self%vcd < 0) self%vcd=0

   end subroutine ReadSatObs
!-------------------------------------------------------------------------

   subroutine SatObsDone(self)

      implicit none

      class(T_SatObs), intent(inout)          ::  self

      if ( allocated(self%clon) )        deallocate(self%clon)
      if ( allocated(self%clat) )        deallocate(self%clat)
      if ( allocated(self%lon) )        deallocate(self%lon)
      if ( allocated(self%lat) )        deallocate(self%lat)
      if ( allocated(self%qa_flag) )     deallocate(self%qa_flag)
      if ( allocated(self%pressure) )    deallocate(self%pressure)
      if ( allocated(self%kernel_trop) ) deallocate(self%kernel_trop)
      if ( allocated(self%amf_trop) )    deallocate(self%amf_trop)
      if ( allocated(self%vcd) )         deallocate(self%vcd)
      if ( allocated(self%vcd_errvar) )  deallocate(self%vcd_errvar)
      if ( allocated(self%nla) )         deallocate(self%nla)
      if ( allocated(self%kernel_trop_m) ) deallocate(self%kernel_trop_m)
      if ( allocated(self%amf_trop_m) )    deallocate(self%amf_trop_m)
      if ( allocated(self%vcd_errvar_m) )  deallocate(self%vcd_errvar_m)
      if ( allocated(self%vcd_m) )         deallocate(self%vcd_m)
      if ( allocated(self%mdlspace_vcd_m) )         deallocate(self%mdlspace_vcd_m)
      if ( allocated(self%mdlspace_vcd) )         deallocate(self%mdlspace_vcd)
      if ( allocated(self%mdlspace_vr) )         deallocate(self%mdlspace_vr)

   end subroutine SatObsDone

!-------------------------------------------------------------------------
end module sat_obs_mod

