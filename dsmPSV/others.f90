subroutine myinput(psvmodel,stationsinf,tlen,imin,imax,r0,r0lat,r0lon,i_source,i_geocentric)
  implicit none
  character(120) :: psvmodel,stationsinf
  real(kind(0d0)) :: tlen,r0,r0lat,r0lon
  integer :: imin,imax,i_source,i_geocentric

  open(unit=1,file='config.slv',status='unknown')
  read(1,*) psvmodel
  read(1,*) stationsinf
  read(1,*) r0,r0lat,r0lon
  ! NF will use rmax instead of 6371 for other planets
  !r0 = 6371.d0-r0 ! because in this version we write the source DEPTH
  read(1,*) tlen
  read(1,*) imin,imax
  read(1,*) i_source
  read(1,*) i_geocentric
  close(1)

  return

end subroutine myinput


subroutine udertorrsgt(icomp,uder,rsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:8)

  ! icomp = 1
  ! This is for P inversion

  if(icomp.eq.1) then ! vertical component
     rsgt(1) = rsgt(1) + uder(1,1)
     rsgt(2) = rsgt(2) + uder(2,1) + uder(1,2)
     rsgt(3) = rsgt(3) + 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(4) = rsgt(4) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  endif

  if(icomp.eq.2) then ! radial component
     rsgt(5) = rsgt(5) + uder(1,1)
     rsgt(6) = rsgt(6) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(7) = rsgt(7) - 5.d-1*uder(2,2) + 5.d-1*uder(3,3)
     rsgt(8) = rsgt(8) - 2.d0*uder(1,2)  - uder(2,1)
  endif
     
  
  return
end subroutine udertorrsgt

!
subroutine udertorsgt(icomp,uder,rsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:2)

  ! icomp = 1
  ! This is for P inversion
  
  rsgt(icomp) = rsgt(icomp) + uder(1,1) + uder(2,2) + uder(3,3)
  
  return
end subroutine udertorsgt

!

subroutine udertotsgt(imt,uder,tsgt)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: uder(1:3,1:3),tsgt(1:4)
  
  ! 1 <= imt <= 4

  ! This is for P inversion

  tsgt(imt) = tsgt(imt) + uder(1,1) + uder(2,2) + uder(3,3)
  return
end subroutine udertotsgt


subroutine udertoStress(uder,stress,A,C,F,L,N)
  implicit none
  complex(kind(0d0)) :: uder(1:3,1:3), stress(1:6)
  real(kind(0d0)) :: A,C,F,L,N
  
  stress(1) = stress(1)+C*uder(1,1)+F*uder(2,2)+F*uder(3,3)
  stress(2) = stress(2)+F*uder(1,1)+A*uder(2,2)+A*uder(3,3)-2.d0*N*uder(3,3)
  stress(3) = stress(3)+F*uder(1,1)+A*uder(2,2)+A*uder(3,3)-2.d0*N*uder(2,2)
  stress(4) = stress(4)+L*(uder(1,2)+uder(2,1))
  stress(5) = stress(5)+L*(uder(1,3)+uder(3,1))
  stress(6) = stress(6)+N*(uder(2,3)+uder(3,2))
  
  return

end subroutine udertoStress
  

!

subroutine setmt(imt,mt)
  implicit none
  integer, intent(in) :: imt
  real(kind(0d0)), dimension(3,3), intent(out) :: mt

  ! We use the  rr, tt, pp, rt, rp, tp order here!!  
  mt = 0.d0
  if (imt==1) mt(1,1) = 1.d0
  if (imt==2) mt(2,2) = 1.d0
  if (imt==3) mt(3,3) = 1.d0
  if (imt==4) mt(1,2) = 1.d0
  if (imt==5) mt(1,3) = 1.d0
  if (imt==6) mt(2,3) = 1.d0

  return
end subroutine setmt


subroutine setf(iif,f)
  implicit none
  integer, intent(in) :: iif
  real(kind(0d0)), dimension(3), intent(out) :: f

  ! We use the  rr, tt, pp, rt, rp, tp order here!!  
  f = 0.d0
  if (iif==1) f(1) = 1.d0
  if (iif==2) f(2) = 1.d0
  if (iif==3) f(3) = 1.d0

  return
end subroutine setf

!

subroutine setmttest(imt,mt)
  implicit none
  ! We use the  rr, tt, pp, rt, rp, tp order here!!  

  real(kind(0d0)) :: mt(3,3)
  integer :: imt
  mt = 0.d0
  if(imt.eq.1) then
     mt(1,1) = 1.d0
     mt(2,2) = 1.d0
     mt(3,3) = 1.d0
  endif

  return
end subroutine setmttest


!

subroutine locallyCartesianDerivatives (u,udr,udt,udp,uder,r,theta)
  implicit none
  complex(kind(0d0)):: u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)
  real(kind(0d0)) :: r,theta
  real(kind(0d0)) :: thetasin,thetacot

  thetasin = sin(theta)
  thetacot = cos(theta)/thetasin

  ! 1,2,3: r,theta,phi; , denotes the partial derivatives

  uder(1,1) = udr(1)
  uder(1,2) = (udt(1)-u(2))/cmplx(r)
  uder(1,3) = (udp(1)/cmplx(thetasin)-u(3))/cmplx(r)

  uder(2,1) = udr(2)
  uder(2,2) = (udt(2)+u(1))/cmplx(r)
  uder(2,3) = (udp(2)/cmplx(thetasin)-u(3)*cmplx(thetacot))/cmplx(r)

  uder(3,1) = udr(3)
  uder(3,2) = udt(3)/cmplx(r)
  uder(3,3) = (udp(3)/cmplx(thetasin)+u(1)+u(2)*cmplx(thetacot))/cmplx(r)

  return
end subroutine locallyCartesianDerivatives

!

subroutine calnl(nzone,vs,iphase,nsl,nll)
  ! counting of nsl and nll.
  implicit none
  integer, intent(in) :: nzone
  real(kind(0d0)), dimension(4,nzone), intent(in) :: vs
  integer, dimension(nzone), intent(out) :: iphase
  integer, intent(out) :: nsl,nll
  integer :: i

  nsl = 0
  nll = 0
  do i = 1,nzone
    if ((vs(1,i)==0.d0).and.(vs(2,i)==0.d0).and.(vs(3,i)==0.d0).and.(vs(4,i)==0.d0)) then
      nll = nll+1
      iphase(i) = 2
    else
      nsl = nsl+1
      iphase(i) = 1
    endif
  enddo

  return
end subroutine calnl

subroutine calgrid(nzone,vrmin,vrmax,vp,vs,rmin,rmax,imax,lmin,tlen,vmin,gridpar,dzpar)
  implicit none
  integer, intent(in) :: nzone,imax,lmin
  real(kind(0d0)), intent(in) :: rmin,rmax,tlen
  real(kind(0d0)), dimension(nzone), intent(in) :: vrmin,vrmax
  real(kind(0d0)), dimension(4,nzone), intent(in) :: vp,vs
  real(kind(0d0)), dimension(nzone), intent(out) :: vmin,gridpar,dzpar
  integer:: i,j,izone
  real(kind(0d0)):: coef1,coef2,vs1,vs2,rh,omega,amax,gtmp
  real(kind(0d0)), dimension(4) :: v
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0

  do izone = 1,nzone
    !!! computing the S-velocity at each zone
    if (vs(1,izone)==0.d0) then
      do i = 1,4
        v(i) = vp(i,izone)
      enddo
    else
      do i = 1,4
        v(i) = vs(i,izone)
      enddo
    endif
    vs1 = 0.d0
    vs2 = 0.d0
    do j = 1,4
      if (j==1) then
        coef1 = 1.d0
      else
        coef1 = coef1*(vrmin(izone)/rmax)
      endif
      if (j==1)then
        coef2 = 1.d0
      else
        coef2 = coef2*(vrmax(izone)/rmax)
      endif
      vs1 = vs1+v(j)*coef1
      vs2 = vs2+v(j)*coef2
    enddo
    rh = vrmax(izone)-vrmin(izone)
    ! computing omega,amax
    omega = 2.d0*pi*dble(imax)/tlen
    if (vs1>=vs2) then
      vmin(izone) = vs2
    else
      vmin(izone) = vs1
    endif
    amax = vrmax(izone)
    gtmp = (omega*omega)/(vmin(izone)*vmin(izone))-((dble(lmin)+0.5d0)*(dble(lmin)+0.5d0))/ &
      (amax*amax)
    if (gtmp>0.d0) then
      dzpar(izone) = dsqrt(1.d0/gtmp)
      gridpar(izone) = rh/dzpar(izone)
    else
      dzpar(izone) = 0.d0
      gridpar(izone) = 0.d0
    endif
  enddo
  ! rearangement of gridpar
  gtmp = 0.d0
  do izone = 1,nzone
    gtmp = gtmp+gridpar(izone)
  enddo
  do izone = 1,nzone
    if (gridpar(izone)>0.d0) then
      gridpar(izone) = gridpar(izone)/gtmp
    else
      rh = vrmax(izone)-vrmin(izone)
      gridpar(izone) = rh/(rmax-rmin)*0.1d0
    endif
  enddo
! re-rearangement of gridpar
  gtmp = 0.d0
  do izone = 1,nzone
    gtmp = gtmp+gridpar(izone)
  enddo
  do izone = 1,nzone
    gridpar(izone) = gridpar(izone)/gtmp
  enddo

  return
end subroutine calgrid


subroutine calra( maxnlay,maxnslay,maxnllay,maxnzone,maxnstack,nlayer,inlayer,jnlayer,jnslay &
  ,jnllay,gridpar,dzpar,nzone,vrmin,vrmax,iphase,rmin,rmax,nslay,nllay,nnl,ra,re, nsta,rsta &
  ,rrsta,istazone,iista,r0,cista) 
  ! Computing the number and the location of grid points.
  
  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0 
  integer:: maxnlay,maxnslay,maxnllay,maxnzone,maxnstack
  integer:: nlayer,inlayer,jnlayer,jnslay,jnllay
  integer:: nzone,iphase(*),nslay,nllay,nnl(maxnzone)
  real(kind(0d0)):: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,r0
  real(kind(0d0)):: ra(maxnlay+maxnzone+1)
  integer:: izone,itmp,i,ntmp
  real(kind(0d0)):: rh,re    
  integer:: nsta
  real(kind(0d0)):: rsta(maxnstack),rrsta(3,maxnstack)
  real(kind(0d0)):: ctmp               ! distance betwee source and the nearst
  integer:: istazone(maxnstack)
  integer:: iista(3,maxnstack)
  integer:: ista,j,cista

  ctmp = 7000.d0
  ! Initializing the data
  nslay = 0
  nllay = 0
  inlayer = 0
  do i=1,maxnlay+maxnzone+1
     ra(i) = 0.d0
  enddo
  do izone=1,nzone
     nnl(izone) = 0
  enddo
  jnlayer = 0
  jnslay = 0
  jnllay = 0
  do i=1,maxnstack
     do j=1,3
        rrsta(j,i) = 0.d0
        iista(j,i) = 0
     enddo
  enddo
  
  do i=1,maxnstack
     istazone(i) = 0
  enddo
  ! computing the number and the location of the grid points
  ra(1) = rmin
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     if(dzpar(izone).eq.0.d0) then
        ntmp = 1
     else
        ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) / 2.d0 / pi  / 7.d-1 + 1 )
     endif
     !                            ! ntmp (see Geller & Takeuchi 1995 6.2)
     nnl(izone) = ntmp
     if ( nnl(izone).lt.5 ) nnl(izone)=5
     if ( iphase(izone).eq.1 ) nslay = nslay + nnl(izone)
     if ( nslay.gt.maxnslay )  stop  'nslay is too large. (calra)'
     if ( iphase(izone).eq.2 ) nllay = nllay + nnl(izone)
     if ( nllay.gt.maxnllay )  stop  'nllay is too large. (calra)'
     do I=1,nnl(izone)
        itmp = itmp + 1
        if ( itmp.gt.maxnlay ) stop  'nlay is too large. (calra)'
        ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo
  
  itmp = 1
  do izone=1,nzone
     do i=1,nnl(izone)
        do ista=1,nsta
           if( (ra(itmp).lt.rsta(ista)).and.(rsta(ista).le.ra(itmp+1)) ) then
              if(i.ne.nnl(izone)) then
                 istazone(ista) = izone
                 if(iphase(istazone(ista)).eq.2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)
                     
                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 if(iphase(istazone(ista)).eq.2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp-1) 
                 rrsta(2,ista) = ra(itmp) 
                 rrsta(3,ista) = ra(itmp+1)
                     
                 iista(1,ista) = i - 1 
                 iista(2,ista) = i 
                 iista(3,ista) = i + 1
              endif
              if(dabs(r0-rsta(ista)).lt.ctmp) then
                 cista = ista
                 ctmp = dabs(r0-rsta(ista))
              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo

  ! recouting the total number of grid points
  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  jnlayer = jnlayer + inlayer
  jnslay  = jnslay  + nslay
  jnllay  = jnllay  + nllay
  
  return
end subroutine calra

!


subroutine calra2(maxnlay,maxnzone,maxnstack,nlayer,inlayer,jnlayer,jnslay,jnllay,gridpar,nzone &
  ,vrmin,vrmax,iphase,rmin,rmax,r0,nslay,nllay,nnl,ra, nsta,rsta,rrsta,istazone,iista) 
  ! Computing the number and the location of grid points.
  
  implicit none
  integer:: maxnlay,maxnzone,maxnstack
  integer:: nlayer,inlayer,jnlayer,jnslay,jnllay
  integer:: nzone,iphase(*),nslay,nllay,nnl(maxnzone)
  real(kind(0d0)):: gridpar(*),vrmin(*),vrmax(*),rmin,rmax,r0
  real(kind(0d0)):: ra(maxnlay+maxnzone+1)
  integer:: izone,itmp,i
  real(kind(0d0)):: rh
  
  integer:: nsta
  real(kind(0d0)):: rsta(maxnstack),rrsta(3,maxnstack)
  integer:: istazone(maxnstack)
  integer:: iista(3,maxnstack)
  integer:: ista,j

  ! Initializing the data
  nslay = 0
  nllay = 0
  inlayer = 0
  do i=1,maxnlay+maxnzone+1
     ra(i) = 0.d0
  enddo
  do izone=1,nzone
     nnl(izone) = 0
  enddo
  do i=1,maxnstack
     do j=1,3
        rrsta(j,i) = 0.d0
        iista(j,i) = 0
     enddo
  enddo
  do i=1,maxnstack
     istazone(i) = 0
  enddo
  jnlayer = 0
  jnslay = 0
  jnllay = 0

  !     computing the number and the location of the grid points
  ra(1) = rmin
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     nnl(izone) = dint( dble(nlayer) * gridpar(izone) )+ 1
     if ( nnl(izone).lt.5 ) nnl(izone)=5
     if ( iphase(izone).eq.1 ) nslay = nslay + nnl(izone)
     if ( iphase(izone).eq.2 ) nllay = nllay + nnl(izone)
     do i=1,nnl(izone)
        itmp = itmp + 1
        ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo
  itmp = 1
  do izone=1,nzone
     do i=1,nnl(izone)
        do ista=1,nsta
           if( (ra(itmp).lt.rsta(ista)).and.(rsta(ista).le.ra(itmp+1)) ) then
              if(i.ne.nnl(izone)) then
                 istazone(ista) = izone
                 if(iphase(istazone(ista)).eq.2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)
                     
                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 if(iphase(istazone(ista)).eq.2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp-1) 
                 rrsta(2,ista) = ra(itmp) 
                 rrsta(3,ista) = ra(itmp+1)
                      
                 iista(1,ista) = i - 1 
                 iista(2,ista) = i 
                 iista(3,ista) = i + 1
              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo
  ! recouting the total number of grid points
  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  jnlayer = jnlayer + inlayer
  jnslay  = jnslay  + nslay
  jnllay  = jnllay  + nllay
  
  return
end subroutine calra2

!

subroutine calsp(nzone,ndc,nsl,nll,iphase,nlayer,nllay,isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr)
  implicit none
  integer, intent(in) :: nzone,ndc,nsl,nll,nllay
  integer, dimension(nzone), intent(in) :: nlayer,iphase
  integer, intent(out) :: isdr,jsdr,ildr,jdr,kdr
  integer, dimension(nzone), intent(out) :: isp,jsp,ksp,issp,ilsp,lsp,jssp
  integer :: i,isl,ill
  ! Initialization of the data
  do i = 1,nzone
    isp(i) = 0
    jsp(i) = 0
    ksp(i) = 0
    issp(i) = 0
    ilsp(i) = 0
    lsp(i) = 0
    jssp(i) = 0
  enddo
  ! computation of isp,jsp,ksp,issp,ilsp,lsp
  isp(1) = 1
  jsp(1) = 1
  ksp(1) = 1
  issp(1) = 1
  ilsp(1) = 1
  lsp(1) = 1
  jssp(1) = 1
  isl = 0
  ill = 0
  do i = 1,ndc
    isp(i+1) = isp(i)+nlayer(i)
    if (iphase(i)==1) then
      jsp(i+1) = jsp(i)+16*nlayer(i)
      ksp(i+1) = ksp(i)+2*(nlayer(i)+1)
      lsp(i+1) = lsp(i)+4*nlayer(i)
      isl = isl+1
      if (isl/=nsl) then
        issp(isl+1) = issp(isl)+4*nlayer(i)
        jssp(isl+1) = jssp(isl)+nlayer(i)+1
      endif
    else
      jsp(i+1) = jsp(i)+4*nlayer(i)
      ksp(i+1) = ksp(i)+(nlayer(i)+1)
      lsp(i+1) = lsp(i)+2*nlayer(i)
      ill = ill+1
      if (ill/=nll) ilsp(ill+1) = ilsp(ill)+4*nlayer(i)
    endif
  enddo
  isdr = issp(nsl)-1+4*nlayer(ndc+1)
  jsdr = jssp(nsl)-1+nlayer(ndc+1)+1
  ildr = 4*nllay
  jdr = jsp(ndc+1)-1+16*nlayer(ndc+1)
  kdr = ksp(ndc+1)-1+2*(nlayer(ndc+1)+1)

  return
end subroutine calsp


subroutine calspo(nlay,nzone,vrmax,iphase,inlayer,ra,rmin,rmax,r0,isp,spo,spn)
  !!! computing the source location.
  implicit none
  integer, intent(in) :: nlay,nzone,inlayer
  real(kind(0d0)), dimension(nzone), intent(in) :: vrmax
  integer, dimension(nzone), intent(in) :: iphase,isp
  real(kind(0d0)), intent(in) :: rmin,rmax
  real(kind(0d0)), dimension(nlay+nzone+1), intent(in) :: ra
  real(kind(0d0)), intent(inout) :: r0
  integer, intent(out) :: spn
  real(kind(0d0)), intent(out) :: spo
  integer :: itmp

  !!! checking the parameter
  if ((r0<rmin).or.(r0>rmax)) stop "The source location is improper.(calspo)"
  spo = 0.d0
  !!! computing spo
  if (r0==rmax) then
    spo = dble(inlayer)-0.01d0
    r0 = ra(inlayer)+(spo-dble(inlayer-1))*(ra(inlayer+1)-ra(inlayer))
  else
    itmp = 2
    do while (r0>=ra(itmp))
      itmp = itmp+1
    enddo
    spo = dble(itmp-2)+(r0-ra(itmp-1))/(ra(itmp)-ra(itmp-1))
    !!! temporal handling
    if ((spo-dble(itmp-2))<0.01d0) then
      spo = dble(itmp-2)+0.01d0
      r0 = ra(itmp-1)+(spo-dble(itmp-2))*(ra(itmp)-ra(itmp-1))
    endif
    if ((spo-dble(itmp-2))>0.99d0) then
      spo = dble(itmp-2)+0.99d0
      r0 = ra(itmp-1)+(spo-dble(itmp-2))*(ra(itmp)-ra(itmp-1))
    endif
  endif
  !!! computing spn
  itmp = 1
  do while (r0>vrmax(itmp)) 
    itmp = itmp+1
  enddo
  spn = itmp
  if (iphase(itmp)/=1) stop "The source is in the liquid layer.(calspo)"
  !!! changing spo
  spo = spo-dble(isp(spn)-1)

  return
end subroutine calspo


subroutine calstg(nlay,nzone,spn,rmax,r0,ra,nnl,rrho,vpv,vph,vsv,vsh,eta,vra,rho,kappa,mu,ecKx &
  ,ecKy,ecKz,ecL,ecN,ecC0,ecF0,ecL0,vnp)
  ! Computing the structure grid points.
  implicit none
  integer, intent(in) :: nlay,nzone,spn
  real(kind(0d0)), dimension(4,nzone), intent(in) :: rrho,vpv,vph,vsv,vsh,eta
  integer, dimension(nzone), intent(in) :: nnl
  real(kind(0d0)), dimension(nlay+nzone+1), intent(in) :: ra
  real(kind(0d0)), intent(in) :: rmax,r0
  real(kind(0d0)), dimension(nlay+2*nzone+1), intent(out) :: vra,rho,kappa,mu
  real(kind(0d0)), dimension(nlay+2*nzone+1), intent(out) :: ecKx,ecKy,ecKz,ecL,ecN
  integer, intent(out) :: vnp
  real(kind(0d0)), intent(out) :: ecC0,ecF0,ecL0
  integer:: izone,i,j,itmp,jtmp
  real(kind(0d0)) :: ecA,ecC,ecF,ecA0
  real(kind(0d0)) :: trho,tvpv,tvph,tvsv,tvsh,teta,coef

  ! initializing the data
  call vecinit(nlay+2*nzone+1,vra)
  call vecinit(nlay+2*nzone+1,rho)
  call vecinit(nlay+2*nzone+1,kappa)
  call vecinit(nlay+2*nzone+1,ecKx)
  call vecinit(nlay+2*nzone+1,ecKy)
  call vecinit(nlay+2*nzone+1,ecKz)
  call vecinit(nlay+2*nzone+1,mu)
  call vecinit(nlay+2*nzone+1,ecL)
  call vecinit(nlay+2*nzone+1,ecN)
  ! computing the structure grid points
  itmp = 0
  jtmp = 0
  do izone = 1,nzone
    do i = 1,nnl(izone)+1
      itmp = itmp+1
      jtmp = jtmp+1
      vra(itmp) = ra(jtmp)
      ! --- evaluating the density and elastic constants at this point
      trho = 0.d0
      tvpv = 0.d0
      tvph = 0.d0
      tvsv = 0.d0
      tvsh = 0.d0
      teta = 0.d0
      do j = 1,4
        if (j==1) then
          coef = 1.d0
        else
          coef = coef*(vra(itmp)/rmax)
        endif
        trho  = trho+rrho(j,izone)*coef
        tvpv  = tvpv+vpv(j,izone)*coef
        tvph  = tvph+vph(j,izone)*coef
        tvsv  = tvsv+vsv(j,izone)*coef
        tvsh  = tvsh+vsh(j,izone)*coef
        teta  = teta+eta(j,izone)*coef
      enddo
      rho(itmp) = trho
      ecL(itmp) = rho(itmp)*tvsv*tvsv
      ecN(itmp) = rho(itmp)*tvsh*tvsh
      ecA = trho*tvph*tvph
      ecC = trho*tvpv*tvpv
      ecF = teta*(ecA-2.d0*ecL(itmp))
      kappa(itmp) = (4.d0*ecA+ecC+4.d0*ecF-4.d0*ecN(itmp))/9.d0
      ecKx(itmp) = ecA-4.d0/3.d0*ecN(itmp)
      ecKy(itmp) = ecF+2.d0/3.d0*ecN(itmp)
      ecKz(itmp) = (ecC+2.d0*ecF)/3.d0
    enddo
    jtmp = jtmp-1
  enddo
  vnp = itmp
  trho = 0.d0
  tvpv = 0.d0
  tvph = 0.d0
  tvsv = 0.d0
  tvsh = 0.d0
  teta = 0.d0
  do j = 1,4
    if (j==1) then
      coef = 1.d0
    else
      coef = coef*(r0/rmax)
    endif
    trho  = trho+rrho(j,spn)*coef
    tvpv  = tvpv+vpv(j,spn)*coef
    tvph  = tvph+vph(j,spn)*coef
    tvsv  = tvsv+vsv(j,spn)*coef
    tvsh  = tvsh+vsh(j,spn)*coef
    teta  = teta+eta(j,spn)*coef
  enddo
  ecL0 = trho*tvsv*tvsv
  ecA0 = trho*tvph*tvph
  ecC0 = trho*tvpv*tvpv
  ecF0 = teta*(ecA0-2.d0*ecL0)

  return
end subroutine calstg


subroutine caltstg(nlay,nzone,rrho,vpv,vph,vsv,vsh,eta,nnl,ra,rmax,tvra,tkappa,tecKx,tecKy, &
  tecKz,tmu,tecL,tecN)
  ! Computing the structure grid points.
  implicit none
  integer, intent(in) :: nlay,nzone
  real(kind(0d0)), dimension(4,nzone), intent(in) :: rrho,vpv,vph,vsv,vsh,eta
  integer, dimension(nzone), intent(in) :: nnl
  real(kind(0d0)), dimension(nlay+nzone+1), intent(in) :: ra
  real(kind(0d0)), intent(in) :: rmax
  real(kind(0d0)), dimension(nlay+2*nzone+1), intent(out) :: tvra,tkappa,tmu
  real(kind(0d0)), dimension(nlay+2*nzone+1), intent(out) :: tecKx,tecKy,tecKz,tecL,tecN
  integer:: izone,i,j,itmp,jtmp
  real(kind(0d0)) :: ecA,ecC,ecF
  real(kind(0d0)) :: trho,tvpv,tvph,tvsv,tvsh,teta,coef

  call vecinit(nlay+2*nzone+1,tvra)
  call vecinit(nlay+2*nzone+1,tkappa)
  call vecinit(nlay+2*nzone+1,tecKx)
  call vecinit(nlay+2*nzone+1,tecKy)
  call vecinit(nlay+2*nzone+1,tecKz)
  call vecinit(nlay+2*nzone+1,tmu)
  call vecinit(nlay+2*nzone+1,tecL)
  call vecinit(nlay+2*nzone+1,tecN)
  !!! computing the structure grid points
  itmp = 0
  jtmp = 0
  do izone = 1,nzone
    do i = 1,nnl(izone)+1
      itmp = itmp+1
      jtmp = jtmp+1
      tvra(itmp) = ra(jtmp)
      !!! evaluating the density and elastic constants at this point
      trho = 0.d0
      tvpv = 0.d0
      tvph = 0.d0
      tvsv = 0.d0
      tvsh = 0.d0
      teta = 0.d0
      do j = 1,4
        if (j==1) then
          coef = 1.d0
        else
          coef = coef*(tvra(itmp)/rmax)
        endif
        trho = trho+rrho(j,izone)*coef
        tvpv = tvpv+vpv(j,izone)*coef
        tvph = tvph+vph(j,izone)*coef
        tvsv = tvsv+vsv(j,izone)*coef
        tvsh = tvsh+vsh(j,izone)*coef
        teta = teta+eta(j,izone)*coef
      enddo
      tecL(itmp) = trho*tvsv*tvsv
      tecN(itmp) = trho*tvsh*tvsh
      ecA = trho*tvph*tvph
      ecC = trho*tvpv*tvpv
      ecF = teta*(ecA-2.d0*tecL(itmp))
      tkappa(itmp) = (4.d0*ecA+ecC+4.d0*ecF-4.d0*tecN(itmp))/9.d0
      tecKx(itmp) = ecA-4.d0/3.d0*tecN(itmp)
      tecKy(itmp) = ecF+2.d0/3.d0*tecN(itmp)
      tecKz(itmp) = (ecC+2.d0*ecF)/3.d0
    enddo
    jtmp = jtmp-1
  enddo

  return
end subroutine caltstg


subroutine calinv(nlay,nzone,vnp,rho,kappa,rhoinv,kappainv)
!!! Computing the inverse of density and elastic constant
  implicit none
  integer, intent(in) :: nlay,nzone,vnp
  real(kind(0d0)), dimension(nlay+2*nzone+1), intent(in) :: rho,kappa
  real(kind(0d0)), dimension(nlay+2*nzone+1), intent(out) :: rhoinv,kappainv
  integer :: i

  do i = 1,vnp
    rhoinv(i) = 1.d0/rho(i)
    kappainv(i) = 1.d0/kappa(i)
  enddo

  return
end subroutine calinv



subroutine submat( nlayer,ha,hb,h )

  ! Subtracting matrix `hb' from matrix `ha'.
  implicit none
  integer:: nlayer
  real(kind(0d0)):: ha(*),hb(*),h(*)
  integer:: i
  
  do i=1,4*nlayer
     h(i) = ha(i) - hb(i)
  enddo
  return
end subroutine submat


!


subroutine calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )

  implicit none
  integer:: maxnzone,nzone,iphase(*)
  integer:: nlayer(maxnzone),jjdr(*),kkdr(*)
  integer:: izone
     
  jjdr(1) = 1
  kkdr(1) = 1
  do izone=1,nzone-1
     if ( iphase(izone).eq.1 ) then
        jjdr(izone+1) = jjdr(izone) + 16 * nlayer(izone)
        if ( iphase(izone+1).eq.1 ) then
           kkdr(izone+1) = kkdr(izone) + 2 * nlayer(izone)
        else
           kkdr(izone+1) = kkdr(izone) + 2 * ( nlayer(izone)+1 )
        endif
     else
        jjdr(izone+1) = jjdr(izone) + 4 * nlayer(izone)
        if ( iphase(izone+1).eq.1 ) then
           kkdr(izone+1) = kkdr(izone) + ( nlayer(izone)+1 )
        else
           kkdr(izone+1) = kkdr(izone) + nlayer(izone)
        endif
     endif
  enddo

  return
end subroutine calspdr

!



subroutine calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0 
  integer:: l,nzone,sufzone
  real(kind(0d0)):: omega,vrmin(*),vrmax(*),vmin(*),dzpar(*),rmax
  integer:: izone
  real(kind(0d0)):: gtmp,tdzpar
  sufzone = 0
  do izone=1,nzone
     gtmp = (omega*omega)/(vmin(izone)*vmin(izone))-((dble(l)+0.5d0)*(dble(l)+0.5d0))/(vrmax(izone)*vrmax(izone))
     if ( gtmp.gt.0.d0 ) then
        tdzpar = sqrt( 1.d0/gtmp )
     else
        if ( vrmax(izone).gt.rmax*(1-2.d0*pi/(dble(l)+0.50)) )  then
           tdzpar = 0.d0
        else
           sufzone = izone
           tdzpar = 0.d0
        endif
     endif
  enddo

  return
end subroutine calmdr

!


subroutine calu0( c0,bvec,u )
  implicit none
  complex(kind(0d0)):: c0,bvec,u

  u = u + c0 * bvec
  
  return
end subroutine calu0

!
     

subroutine calulcd0( c0,c0der,rsta,theta, bvec,bvecdt,bvecdp,ulcd )
  
  implicit none
  complex(kind(0d0)):: c0,c0der,bvec(3),bvecdt(3),bvecdp(3),ulcd(9)
  real(kind(0d0)):: rsta,theta
  complex(kind(0d0)):: u1,uder11,uder12,uder13
  
  u1 = c0 * bvec(1)
  uder11 = c0der * bvec(1)
  uder12 = c0 * bvecdt(1)
  uder13 = c0 * bvecdp(1)
  
  ulcd(1) = ulcd(1) + uder11
  ulcd(2) = ulcd(2) + uder12 / rsta
  ulcd(3) = ulcd(3) + uder13 / rsta / dsin(theta)
  ulcd(5) = ulcd(5) + u1 / rsta
  ulcd(9) = ulcd(9) + u1 / rsta
   
  return
end subroutine calulcd0

!
     

subroutine calu( c0,lsq,bvec,u )
  implicit none
  real(kind(0d0)):: lsq
  complex(kind(0d0)):: c0(2),bvec(3),u(3)
  
  u(1) = u(1) + c0(1) * bvec(1)
  u(2) = u(2) + c0(2) * bvec(2) / dcmplx(lsq)
  u(3) = u(3) + c0(2) * bvec(3) / dcmplx(lsq)

  return
end subroutine calu


!

subroutine calup(c1,c2,lsq,bvec,u)
  implicit none
  real(kind(0d0)) :: lsq
  complex(kind(0d0)) :: c1,c2, bvec(1:3), u(1:3)

  u(1) = u(1) + c1*bvec(1)
  u(2) = u(2) + c2*bvec(2)/dcmplx(lsq)
  u(3) = u(3) + c2*bvec(3)/dcmplx(lsq)
  
  return
end subroutine calup

!

subroutine calup0(c1,bvec,u)
  implicit none
  complex(kind(0d0)) :: c1,u(1:3),bvec(1:3)

  u(1) = u(1) + c1*bvec(1)
  
  return
end subroutine calup0


!


subroutine calulcd( c0,c0der,lsq,rsta,theta, bvec,bvecdt,bvecdp,ulcd )

  implicit none
  real(kind(0d0)):: lsq,rsta,theta
  complex(kind(0d0)):: c0(2),c0der(2)
  complex(kind(0d0)):: bvec(3),bvecdt(3),bvecdp(3),ulcd(9)
  
  complex(kind(0d0)):: u1,u2,u3
  complex(kind(0d0)):: uder11,uder12,uder13
  complex(kind(0d0)):: uder21,uder22,uder23
  complex(kind(0d0)):: uder31,uder32,uder33

  u1 = c0(1) * bvec(1)
  u2 = c0(2) * bvec(2) / dcmplx(lsq)
  u3 = c0(2) * bvec(3) / dcmplx(lsq)
  ! partial derivatives of u
  uder11 = c0der(1) * bvec(1)
  uder12 = c0(1) * bvecdt(1)
  uder13 = c0(1) * bvecdp(1)
  uder21 = c0der(2) * bvec(2) / dcmplx(lsq)
  uder22 = c0(2) * bvecdt(2) / dcmplx(lsq)
  uder23 = c0(2) * bvecdp(2) / dcmplx(lsq)
  uder31 = c0der(2) * bvec(3) / dcmplx(lsq)
  uder32 = c0(2) * bvecdt(3) / dcmplx(lsq)
  uder33 = c0(2) * bvecdp(3) / dcmplx(lsq)
  ! locally Cartesian derivatives of u
  ulcd(1) = ulcd(1) + uder11
  ulcd(2) = ulcd(2) + ( uder12 - u2 ) / rsta 
  ulcd(3) = ulcd(3) + ( uder13 / dsin(theta) - u3 ) / rsta
  ulcd(4) = ulcd(4) + uder21
  ulcd(5) = ulcd(5) + ( uder22 + u1 ) / rsta
  ulcd(6) = ulcd(6) + ( uder23 - u3 * dcos(theta) )  / rsta / dsin(theta)
  ulcd(7) = ulcd(7) + uder31
  ulcd(8) = ulcd(8) + uder32 / rsta
  ulcd(9) = ulcd(9)  + ( ( uder33 + u2 * dcos(theta) ) / dsin(theta) + u1 ) / rsta
  return
end subroutine calulcd
    
!


subroutine matinit( n1,n2,a )
  implicit none
  integer:: n1,n2,i,j
  real(kind(0d0)):: a(n1,*)

  do j=1,n2
     do i=1,n1
        a(i,j) = 0.d0
     enddo
  enddo
  
  return
end subroutine matinit

!


subroutine cmatinit( n1,n2,a )
  implicit none
  integer:: n1,n2,i,j
  complex(kind(0d0)):: a(n1,*)
  
  do j=1,n2
     do i=1,n1
        a(i,j) = dcmplx( 0.d0 )
     enddo
  enddo
  return
end subroutine cmatinit

!
     

subroutine vecinit( nn,b )

  ! Filling zero to the vector 'g'.
  implicit none
  integer:: nn,i
  real(kind(0d0)):: b(*)
     
  do i=1,nn
     b(i) = 0.d0
  enddo
  return
end subroutine vecinit

!


subroutine cvecinit( nn,b )

  ! Filling zero to the vector 'g'.
  implicit none 
  integer:: nn,i
  complex(kind(0d0)):: b(*)
  
  do i=1,nn
     b(i) = dcmplx( 0.d0 )
  enddo
  return
end subroutine cvecinit

!


subroutine interpolate( ncomp,nderiv,rsta,rrsta,g,u )

  implicit none
  integer:: ncomp,nderiv
  real(kind(0d0)):: rsta,rrsta(3)
  complex(kind(0d0)):: g(3*ncomp),u(ncomp)
  real(kind(0d0)):: dh(3)
  
  integer:: ip(3),ier,i,itmp,icomp
  complex(kind(0d0)):: a(3,3),b(3),wk(3)
  real(kind(0d0)):: eps
  eps = -1.d0
  
  do icomp=1,ncomp
     u(icomp) = dcmplx(0.d0)
  enddo
  do i=1,3
     dh(i) = rrsta(i) - rsta
  enddo
   
  if( (dh(2).eq.0.d0).and.(nderiv.eq.0)) then
     itmp = ncomp + 1
     do icomp=1,ncomp
        u(icomp) = g(itmp)
        itmp = itmp + 1
     enddo
     return
  endif    
  do i=1,3
     a(1,i) = dcmplx( 1.d0 )
     a(2,i) = dcmplx( dh(i) )
     a(3,i) = dcmplx( dh(i) * dh(i) / 2.d0 )
  enddo
  call fillinpb(nderiv,b)
  call glu(a,3,3,b,eps,wk,ip,ier)
    
  
  do icomp=1,ncomp
     do i=1,3
        u(icomp) = u(icomp) + b(i) * g( ncomp * (i-1) + icomp )
     enddo
  enddo
end subroutine interpolate

!


subroutine fillinpb( nderiv,b )

  implicit none
  integer:: nderiv
  complex(kind(0d0)):: b(3)
  
  if( (nderiv.ne.0).and.(nderiv.ne.1).and.(nderiv.ne.2) ) stop 'invalid argument (fillinpb)'
  if(nderiv.eq.0) then
     b(1) = dcmplx( 1.d0 )
     b(2) = dcmplx( 0.d0 )
     b(3) = dcmplx( 0.d0 )
  elseif(nderiv.eq.1) then
     b(1) = dcmplx( 0.d0 )
     b(2) = dcmplx( 1.d0 )
     b(3) = dcmplx( 0.d0 )
  elseif(nderiv.eq.2) then
     b(1) = dcmplx( 0.d0 )
     b(2) = dcmplx( 0.d0 )
     b(3) = dcmplx( 1.d0 )
  endif
  
  return
end subroutine fillinpb


!



subroutine calamp( g,l,lsuf,maxamp,ismall,ratl )
  implicit none
  integer:: ismall,l,lsuf
  real(kind(0d0)):: maxamp,ratl
  complex(kind(0d0)):: g(2)
  real(kind(0d0)):: amp,ampratio
  
  ampratio = 0.d0
  amp = dsqrt( zabs( g(1) )**2 + zabs( g(2) )**2 )
  if ( amp.gt.maxamp ) maxamp = amp
  if ( (amp.ne.0.d0).and.(maxamp.ne.0.d0) ) ampratio = amp / maxamp
  if ( ( ampratio.lt.ratl ).and.( l.ge.lsuf ) ) then
     ismall = ismall + 1
  else
     ismall = 0
  endif

  return
end subroutine calamp

!


subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
  implicit none
  integer:: nzone,lsuf
  real(kind(0d0)):: omega,vrmax(*),vsv(4,*)
  real(kind(0d0)):: tvs,coef
  integer:: i

  tvs = 0.d0
  do i=1,4
     if(i.eq.1) then
        coef = 1.d0
     else
        coef = coef 
     endif
     tvs = tvs + ( vsv(i,nzone) ) * coef
  enddo
  lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
  return
end subroutine callsuf

subroutine calra_psv(nzone,vrmin,vrmax,iphase,dzpar,re,nlayer,nslay,nllay,nnl)
  implicit none
  integer, intent(in) :: nzone
  real(kind(0d0)), intent(in) :: re
  integer, dimension(nzone), intent(in) :: iphase
  real(kind(0d0)), dimension(nzone), intent(in) :: vrmin,vrmax,dzpar
  integer, intent(out) :: nlayer,nslay,nllay
  integer, dimension(nzone), intent(out) :: nnl
  integer :: izone
  real(kind(0d0)) :: rh
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0

  nlayer = 0
  nslay = 0
  nllay = 0
  nnl = 0
  do izone = 1,nzone
    rh = vrmax(izone)-vrmin(izone)
    !!! see Geller & Takeuchi 1995 6.2
    if (dzpar(izone)/=0.d0) nnl(izone) = int(sqrt(3.3d0/re)*rh/dzpar(izone)/2.d0/pi/7.d-1+1)
    if (nnl(izone)<5) nnl(izone) = 5
    nlayer = nlayer+nnl(izone)
    if (iphase(izone)==1) nslay = nslay+nnl(izone)
    if (iphase(izone)==2) nllay = nllay+nnl(izone)
  enddo

  return
end subroutine calra_psv


subroutine calra2_psv(nlayer,nzone,nsta,rmin,rs,nnl,iphase,vrmin,vrmax,rsta,cista,ciista,istazone,iista,ra,rrsta)
  implicit none
  integer, intent(in) :: nlayer,nzone,nsta
  real(kind(0d0)), intent(in) :: rmin,rs
  integer, dimension(nzone), intent(in) :: nnl,iphase
  real(kind(0d0)), dimension(nzone), intent(in) :: vrmin,vrmax
  real(kind(0d0)), dimension(nsta), intent(in) :: rsta
  integer, intent(out) :: cista,ciista
  integer, dimension(nsta), intent(out) :: istazone
  integer, dimension(3,nsta), intent(out) :: iista
  real(kind(0d0)), dimension(nlayer+nzone+1), intent(out) :: ra
  real(kind(0d0)), dimension(3,nsta), intent(out) :: rrsta
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: izone,itmp,i,ista
  real(kind(0d0)) :: rh,ctmp

  ra = 0
  ra(1) = rmin
  ciista = 0
  ctmp = 6371.d0
  itmp = 1
  do izone = 1,nzone
    rh = vrmax(izone)-vrmin(izone)
    do i = 1,nnl(izone)
      itmp = itmp+1
      ra(itmp) = vrmin(izone)+rh*dble(i)/dble(nnl(izone))
    enddo
  enddo
  itmp = 1
  do izone = 1,nzone
    do i = 1,nnl(izone)
      do ista = 1,nsta
        if ((ra(itmp)<rsta(ista)).and.(rsta(ista)<=ra(itmp+1))) then
          if (i/=nnl(izone)) then
            istazone(ista) = izone
            if (iphase(istazone(ista))==2) stop "rsta is in liquid layer (calra2_psv)"
            rrsta(1,ista) = ra(itmp)
            rrsta(2,ista) = ra(itmp+1)
            rrsta(3,ista) = ra(itmp+2)
            iista(1,ista) = i
            iista(2,ista) = i+1
            iista(3,ista) = i+2
          else
            istazone(ista) = izone
            if (iphase(istazone(ista))==2) stop "rsta is in liquid layer (calra2_psv)"
            rrsta(1,ista) = ra(itmp-1)
            rrsta(2,ista) = ra(itmp)
            rrsta(3,ista) = ra(itmp+1)
            iista(1,ista) = i-1
            iista(2,ista) = i
            iista(3,ista) = i+1
          endif
          if ((abs(rs-rsta(ista))<ctmp).and.(abs(rs-rsta(ista))>=0.d0)) then
            cista = ista
            ciista = itmp
            ctmp = abs(rs-rsta(ista))
          endif
        endif
      enddo
      itmp = itmp+1
    enddo
  enddo
  if (cista==0) cista = 1

  return
end subroutine calra2_psv


subroutine calcutd(nzone,nnlayer,nnl,tmpc,rat,nn,iphase,spo,spn, ra,kkdr,kc)
  implicit none
  integer :: nzone,nn,nnlayer,spn,kkdr(1:nzone),kc,iphase(1:nzone),nnl(1:nzone)
  complex(kind(0d0)) :: tmpc(1:nn)
  real(kind(0d0)) :: rat,spo,ra(1:nnlayer+nzone+1)
  integer :: nc
  real(kind(0d0)) :: cU(nn),cV(nn),rc
  real(kind(0d0)) :: maxamp,amp(nn)
  integer :: iz,jz,jj,i,ml(nzone),tzone

  do jj=1,nn
     cU(jj) = 0.d0
     cV(jj) = 0.d0
  enddo
  iz = 2
  jz = 1
  do jj=1,nn
     if(iz.le.nzone) then
        if(jj.eq.kkdr(iz)) then
           if(iphase(iz).ne.iphase(iz-1)) jz = jz - 1
           iz = iz + 1
        endif
     endif
     if(iphase(iz-1).eq.1) then
        if(mod((jj-kkdr(iz-1)),2).eq.1) then ! U
           cU(jz) = cdabs(tmpc(jj))
           jz = jz + 1
        else		! V
        endif
     else ! U in fluid
        cU(jz) = cdabs(tmpc(jj))
        jz = jz + 1
     endif
  enddo

  maxamp = -1.d0
  do i=1,jz-1
     amp(i) = cU(i)
     if(maxamp.lt.amp(i)) maxamp = amp(i)
  enddo
!
  maxamp = maxamp * rat ! threshold value
!
  nc = 1
  do i=1,jz-1
     if(amp(i).gt.maxamp) then
        nc = i
        cycle
     endif
  enddo
  i = 1
  do jj=1,nzone
     i = i + nnl(jj)
     ml(jj) = i
  enddo
  do jj=nzone,1,-1
     if(ml(jj).gt.nc) tzone = jj
  enddo
  rc = ra(nc)
  
  do i=1,jz-1
     if( (ra(i).le.rc).and.(rc.lt.ra(i+1)) ) then
        nc = i
        if(tzone.eq.1) then ! case(tzone is innermost zone)
           if(iphase(tzone).eq.1) kc = 1 + 2 * nc
           if(iphase(tzone).eq.2) kc = 1 + nc
        else 
           if(iphase(tzone).eq.1) then
              kc = kkdr(tzone) + 2 * (nc - ml(tzone-1))
           endif
           if(iphase(tzone).eq.2) then
              kc = kkdr(tzone) + nc - ml(tzone-1)
           endif
        endif
     endif
  enddo
  
  return
end subroutine calcutd




subroutine translat(geodetic,geocentric)

  implicit none
  real(kind(0d0)),parameter ::  flattening = 1.d0 / 298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0 
  real(kind(0d0)) :: geocentric, geodetic 
  real(kind(0d0)) :: tmp
  integer :: flag
  flag = 0
  if(geodetic .gt. 90.d0) then
     geodetic = 1.8d2 - geodetic
     flag = 1
  endif
  
  geodetic = geodetic / 1.8d2 * pi
  geocentric = datan((1.d0-flattening)*(1.d0-flattening)* dtan(geodetic) )
  geocentric = geocentric * 1.8d2 / pi
  
  if(flag .eq. 1) then
     geocentric = 1.8d2 - geocentric
  endif

  return
end subroutine translat


subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
  
  implicit none
  real(kind(0d0)), parameter:: pi = 3.1415926535897932d0 
  
  real(kind(0d0)) ::  ievla,ievlo,istla,istlo
  real(kind(0d0)) ::  evla,evlo,stla,stlo
  real(kind(0d0)) :: theta,phi
  real(kind(0d0)) :: gcarc,az
  real(kind(0d0)) :: tc,ts

  ! transformation to spherical coordinates
  
  evla = 90.d0 - ievla
  stla = 90.d0 - istla
  
  evla = evla / 1.8d2 * pi
  evlo = ievlo / 1.8d2 * pi
  stla = stla / 1.8d2 * pi
  stlo = istlo / 1.8d2 * pi
  
  gcarc = dacos( dcos(evla) * dcos(stla) + dsin(evla) * dsin(stla) * dcos(evlo - stlo) )
  
  tc = (dcos(stla)*dsin(evla)-dsin(stla)*dcos(evla)*dcos(stlo-evlo))/dsin(gcarc)
  ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)

  az = dacos(tc)
  if( ts .lt. 0.d0 ) az = -1.d0 * az
  
  az = az * 1.8d2 / pi
  
  gcarc = gcarc * 1.8d2 / pi

  theta = gcarc
  phi   = 180.d0 - az
  return
end subroutine calthetaphi

subroutine calstg4onedepth(nzone,vrmin,vrmax,rho,vpv,vph,vsv,vsh,eta,rmax,r,updown,ecA,ecC,ecF,ecL,ecN)
  implicit none
  integer, intent(in) :: nzone,updown
  real(kind(0d0)), intent(in) :: rmax,r
  real(kind(0d0)), dimension(nzone), intent(in) :: vrmin,vrmax
  real(kind(0d0)), dimension(4,nzone), intent(in) :: rho,vpv,vph,vsv,vsh,eta
  real(kind(0d0)), intent(out) :: ecA,ecC,ecF,ecL,ecN
  integer :: i,izone,spn
  real(kind(0d0)) :: trho,tvpv,tvph,tvsv,tvsh,teta,coef

  spn = 0
  do izone = 1,nzone
    if ((vrmin(izone)<r).and.(vrmax(izone)>r)) then
      spn = izone
    endif
  enddo
  if (vrmax(nzone)==r) spn = nzone
  if ((vrmin(spn)==r).and.(updown==-1)) then
    spn = spn-1
  endif

  trho = 0.d0
  tvpv = 0.d0
  tvph = 0.d0
  tvsv = 0.d0
  tvsh = 0.d0
  teta = 0.d0
  do i = 1,4
    if (i==1) then
      coef = 1.d0
    else
      coef = coef*(r/rmax)
    endif
    trho = trho+rho(i,spn)*coef
    tvpv = tvpv+vpv(i,spn)*coef
    tvph = tvph+vph(i,spn)*coef
    tvsv = tvsv+vsv(i,spn)*coef
    tvsh = tvsh+vsh(i,spn)*coef
    teta = teta+eta(i,spn)*coef
  enddo
  ecL = trho*tvsv*tvsv
  ecN = trho*tvsh*tvsh
  ecA = trho*tvph*tvph
  ecC = trho*tvpv*tvpv
  ecF = teta*(ecA-2.d0*ecL)

  return
end subroutine calstg4onedepth
