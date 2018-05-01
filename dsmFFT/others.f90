subroutine myinput(nodesfile,normsfile,indexfile,mirrorfile,tlen,imin,imax,fmin,fmax,r0,lat0,lon0,i_source,i_channel,i_geocentric,i_filter,i_basis,mt)
  implicit none
  character(len=40) :: nodesfile,normsfile,indexfile,mirrorfile
  real(kind(0d0)) :: tlen,fmin,fmax,r0,lat0,lon0
  real(kind(0d0)), dimension(6) :: mt
  integer :: imin,imax,i_source,i_channel,i_geocentric,i_filter,i_basis

  open(1,file="config.fft")
  read(1,*) nodesfile
  read(1,*) normsfile
  read(1,*) indexfile
  read(1,*) r0,lat0,lon0
  r0 = 6371.d0-r0
  read(1,*) tlen
  read(1,*) imin,imax
  read(1,*) mt(1),mt(2),mt(3),mt(4),mt(5),mt(6)
  read(1,*) fmin,fmax
  read(1,*) i_source
  read(1,*) i_channel
  read(1,*) i_basis
  read(1,*) i_geocentric
  read(1,*) i_filter
  read(1,*) mirrorfile
  close(1)

  return

end subroutine myinput

subroutine tensorFFT_real(n,imin,np1,cvec,rvec,tlen,i_channel)
  implicit none
  integer :: i,j,n,imin,np1,n1,m1,i_channel
  complex(kind(0d0)), dimension(1:n,0:2*np1-1) :: cvec
  real(kind(0e0)), dimension(1:n,0:2*np1-1) :: rvec
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: omegai,tlen,samplingHz
  complex(kind(0d0)), parameter :: ii = (0.d0, -1.d0)

  omegai = -dlog(1.d-2)/tlen
  samplingHz = dble(2*np1)/tlen
  do j = 1,n
    do i = imin,np1-1
      n1 = np1+i
      m1 = np1-i
      if (i_channel<1) then
        cvec(j,m1) = cvec(j,m1)*ii/(omegai*dble(m1)) !!! velocity to displacement
      elseif (i_channel>1) then
        cvec(j,m1) = cvec(j,m1)/ii*(omegai*dble(m1)) !!! velocity to acceleration
      endif
      cvec(j,n1) = conjg(cvec(j,m1))
    enddo
  enddo
  do j = 1,n
    call cdft(4*np1,cos(pi/(2*np1)),sin(pi/(2*np1)),cvec(j,0:2*np1-1))
    do i = 0,2*np1-1
      cvec(j,i) = dble(cvec(j,i))*dble(exp(omegai*dble(i)/samplingHz))/tlen*1.d3
      rvec(j,i) = real(cvec(j,i))
    enddo
  enddo

  return

end subroutine tensorFFT_real

subroutine lsmoothfinder(tlen,np0,freq,lsmooth)
  implicit none
  real(kind(0d0)) :: tlen,freq
  integer :: np0,np,lsmooth,i

  np = 1
  do while (np<np0)
    np = np*2
  enddo
  lsmooth = int(0.5*tlen*freq/dble(np))
  i = 1
  do while (i<lsmooth)
    i = i*2
  enddo
  lsmooth = i

  return 

end subroutine lsmoothfinder

subroutine bwfilt(x,y,dt,n,irek,norder,f1,f2)
  ! recursive filtering of data with butterworth filter
  ! x: input array
  ! y: output array
  ! dt: time increment
  ! n: number of data points
  ! irek=0: forward filtering only
  ! irek=1: forward and backward filtering
  ! norder: order of butterworth filter
  ! norder=0: only filtering, no determination of coefficients
  ! norder<0: no starplots of transfer function and impulse response
  ! f1: low cutoff frequency (Hz)
  ! f1=0: low pass filter
  ! f2: high cutoff frequency (Hz)
  ! f2>0.5/dt: high pass filter
  implicit none
  real(kind(0d0)), dimension(1) :: x,y
  real(kind(0d0)), dimension(10) :: a,b1,b2
  real(kind(0d0)) :: dt,f1,f2
  integer :: iunit,npoles,norder,irek,n,lx
  
  iunit = 3
  if (norder.ne.0) then
    npoles = iabs(norder)
    !determination of filter coefficients
    call bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    if (norder.ge.0) then
      !plot of transfer function and impuulse response
      lx = 100
      !filtering
    endif
  endif
  if (n.ne.0) then
    call rekurs(x,y,n,a,b1,b2,npoles,irek)
  endif

  return

end subroutine bwfilt

subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag.ne.0: forward and backward filtering
  implicit none
  real(kind(0d0)), dimension(10) :: z,z1,z2,a,b1,b2
  real(kind(0d0)) :: x1,x2
  integer :: ndat,npoles,iflag,n,i
  real(kind(0d0)), dimension(ndat) :: x,y
  
  !forward
  x1 = 0.d0
  x2 = 0.d0
  do i = 1,npoles
    z1(i) = 0.d0
    z2(i) = 0.d0
  enddo
  do n = 1,ndat
    z(1) = a(1)*(x(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
    do i = 2,npoles
      z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
    enddo
    x2 = x1
    x1 = x(n)
    do i = 1,npoles
      z2(i) = z1(i)
      z1(i) = z(i)
    enddo
    y(n) = z(npoles)
  enddo
  if (iflag.eq.0) then
    return
  endif
  !backward
  x1 = 0.d0
  x2 = 0.d0
  do i = 1,npoles
    z1(i) = 0.d0
    z2(i) = 0.d0
  enddo
  do n = ndat,1,-1
    z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
    do i = 2,npoles
      z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
    enddo
    x2 = x1
    x1 = y(n)
    do i = 1,npoles
       z2(i) = z1(i)
       z1(i) = z(i)
    enddo
    y(n) = z(npoles)
  enddo

  return

end subroutine rekurs

subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter
  implicit none
  real(kind(0d0)), dimension(10) :: a,b1,b2
  complex(kind(0d0)), dimension(20) :: s
  complex(kind(0d0)) :: t1,t2,p
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: f1,f2,dt,d2,w0,w1,w2,ssum,sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles
  
  if (npoles.gt.10) stop 'npoles greater than 10: STOP'
  d2 = 2.d0/dt
  w1 = d2*tan(2.d0*pi*f1/d2)
  w2 = d2*tan(2.d0*pi*f2/d2)
  w0 = 0.5*(w2-w1)
  i = 1
  npol2 = npoles/2+1
  do n = 1,npoles
    p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
    t1 = p*cmplx(w0,0.d0)
    t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
    s(i) = t1+t2
    s(i+1) = t1-t2
    i = i+2
  enddo 
  do n = 1,npoles
    ssum = 2*real(s(n))
    sprod = dble(s(n)*conjg(s(n)))
    fact1 = d2*d2-d2*ssum+sprod
    fact2 = 2.d0*(sprod-d2*d2)
    fact3 = d2*d2+d2*ssum+sprod
    a(n) = 2.d0*d2*w0/fact1
    b1(n) = fact2/fact1
    b2(n) = fact3/fact1
  enddo

  return

end subroutine bpcoeff
 
subroutine cdft(n,wr,wi,c)
  implicit none
  integer :: n,i,j,k,l,m
  real(kind(0d0)), dimension(0:n-1) :: a          
  real(kind(0d0)) :: wr,wi,wmr,wmi,wkr,wki 
  real(kind(0d0)) :: wdr,wdi,ss,xr,xi
  complex(kind(0d0)),dimension(0:n/2-1) :: c
  
  do i = 0,n/2-1
    a(2*i) = dble(c(i))
    a(2*i+1) = imag(c(i))
  enddo
  wmr = wr
  wmi = wi
  m = n
  do while (m.gt.4)
    l = m/2
    wkr = 1
    wki = 0
    wdr = 1-2*wmi*wmi
    wdi = 2*wmi*wmr
    ss = 2*wdi
    wmr = wdr
    wmi = wdi
    do j = 0,n-m,m
      i = j+l
      xr = a(j)-a(i)
      xi = a(j+1)-a(i+1)
      a(j) = a(j)+a(i)
      a(j+1) = a(j+1)+a(i+1)
      a(i) = xr
      a(i+1) = xi
      xr = a(j+2)-a(i+2)
      xi = a(j+3)-a(i+3)
      a(j+2) = a(j+2)+a(i+2)
      a(j+3) = a(j+3)+a(i+3)
      a(i+2) = wdr*xr-wdi*xi
      a(i+3) = wdr*xi+wdi*xr
    enddo
    do k = 4,l-4,4
      wkr = wkr-ss*wdi
      wki = wki+ss*wdr
      wdr = wdr-ss*wki
      wdi = wdi+ss*wkr
      do j = k,n-m+k,m
        i = j+l
        xr = a(j)-a(i)
        xi = a(j+1)-a(i+1)
        a(j) = a(j)+a(i)
        a(j+1) = a(j+1)+a(i+1)
        a(i) = wkr*xr-wki*xi
        a(i+1) = wkr*xi+wki*xr
        xr = a(j+2)-a(i+2)
        xi = a(j+3)-a(i+3)
        a(j+2) = a(j+2)+a(i+2)
        a(j+3) = a(j+3)+a(i+3)
        a(i+2) = wdr*xr-wdi*xi
        a(i+3) = wdr*xi+wdi*xr
      enddo
    enddo
    m = l
  enddo
  if (m.gt.2) then
    do j = 0,n-4,4
      xr = a(j)-a(j+2)
      xi = a(j+1)-a(j+3)
      a(j) = a(j)+a(j+2)
      a(j+1) = a(j+1)+a(j+3)
      a(j+2) = xr
      a(j+3) = xi
    enddo
  endif
  if (n.gt.4) call bitrv2(n,a)
  
  do i = 0,n/2-1
    c(i) = dcmplx(a(2*i),a(2*i+1))
  enddo
  
end subroutine cdft

subroutine bitrv2(n,a)
  implicit none
  integer :: n,j,j1,k,k1,l,m,m2,n2
  real(kind(0d0)), dimension(0:n-1) :: a
  real(kind(0d0)) :: xr,xi
  
  m = n/4
  m2 = 2*m
  n2 = n-2
  k = 0
  do j = 0,m2-4,4
    if (j.lt.k) then
      xr = a(j)
      xi = a(j+1)
      a(j) = a(k)
      a(j+1) = a(k+1)
      a(k) = xr
      a(k+1) = xi
    else if (j.gt.k) then
      j1 = n2-j
      k1 = n2-k
      xr = a(j1)
      xi = a(j1+1)
      a(j1) = a(k1)
      a(j1+1) = a(k1+1)
      a(k1) = xr
      a(k1+1) = xi
    endif
    k1 = m2+k
    xr = a(j+2)
    xi = a(j+3)
    a(j+2) = a(k1)
    a(j+3) = a(k1+1)
    a(k1) = xr
    a(k1+1) = xi
    l = m
    do while (k.ge.l)
      k = k-l
      l = l/2
    enddo
    k = k+l
  enddo

  return

end subroutine bitrv2

subroutine translat(geodetic,geocentric)
  implicit none
  real(kind(0d0)), parameter :: flattening = 1.d0/298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0 
  real(kind(0d0)) :: geocentric, geodetic 
  real(kind(0d0)) :: tmp
  integer :: flag

  flag = 0
  if (geodetic.gt.90.d0) then
    geodetic = 1.8d2-geodetic
    flag = 1
  endif
  
  geodetic = geodetic/1.8d2*pi
  geocentric = datan((1.d0-flattening)*(1.d0-flattening)*dtan(geodetic))
  geocentric = geocentric*1.8d2/pi
  
  if (flag.eq.1) then
    geocentric = 1.8d2-geocentric
  endif

  return

end subroutine translat

subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
  implicit none
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0 
  real(kind(0d0)) :: ievla,ievlo,istla,istlo
  real(kind(0d0)) :: evla,evlo,stla,stlo
  real(kind(0d0)) :: theta,phi
  real(kind(0d0)) :: gcarc,az
  real(kind(0d0)) :: tc,ts

  ! transformation to spherical coordinates
  evla = 90.d0-ievla
  stla = 90.d0-istla
  evla = evla/1.8d2*pi
  evlo = ievlo/1.8d2*pi
  stla = stla/1.8d2*pi
  stlo = istlo/1.8d2*pi
  gcarc = dacos(dcos(evla)*dcos(stla)+dsin(evla)*dsin(stla)*dcos(evlo-stlo))
  tc = (dcos(stla)*dsin(evla)-dsin(stla)*dcos(evla)*dcos(stlo-evlo))/dsin(gcarc)
  ts = dsin(stla)*dsin(stlo-evlo)/dsin(gcarc)
  az = dacos(tc)
  if (ts.lt.0.d0) az = -1.d0*az
  az = az*1.8d2/pi
  gcarc = gcarc*1.8d2/pi
  theta = gcarc
  phi = 180.d0-az

  return

end subroutine calthetaphi
