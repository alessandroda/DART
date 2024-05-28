!------------------------------------------------------------
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~ partial columns from local model:
!    y = sum_{l=1,nla} x
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine PartialColumn(nretr,nlayer,nla,x,Sx)
   implicit none
!-in/out ----------------
   integer,intent(in)     :: nlayer
   integer,intent(in)     :: nretr
   integer,intent(in)     :: nla(nretr)
   real,intent(in)        :: x(nlayer)
   real,intent(out)       :: Sx(nretr)
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
   real,intent(in)       :: M(nretr)
   real,intent(in)       :: A(nretr,nlayer)
   real,intent(in)       :: x(nlayer)
   real,intent(in)       :: Sx(nretr)
   real,intent(out)      :: M_m(nretr)
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
!------------------------------------------------------------
