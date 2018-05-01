program  myTraPSV
  implicit none
  character(120) :: psvmodel,stationsinf
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0
  real(kind(0d0)), parameter :: re = 1.d-2
  real(kind(0d0)), parameter :: ratc = 1.d-10
  real(kind(0d0)), parameter :: ratl = 1.d-4
  integer, parameter :: maxlmax = 20000
  integer, dimension(2) :: ltmp 
  real(kind(0d0)) :: tlen,r0,r0lat,r0lon
  real(kind(0d0)), dimension(:), allocatable :: stla,stlo,r_,theta,phi
  real(kind(0d0)), dimension(:), allocatable :: A0sta,C0sta,F0sta,L0sta,N0sta
  integer, dimension(:), allocatable :: updown
  real(kind(0d0)), dimension(:,:), allocatable :: rrsta
  integer, dimension(:,:), allocatable :: iista
  integer :: r_n,ciista,ir_,imt,icomp,idepth,itheta,theta_n,nsta
  character(120) :: coutfile
  integer :: imin,imax
  integer :: i_source,i_geocentric
  ! ---------------------------<< variables >>---------------------------
  ! variable for the trial function
  integer :: nnlayer,nlay
  integer, dimension(:), allocatable :: nlayer
  integer :: nslay,nllay
  integer :: l,m
  real(kind(0d0)), dimension(:), allocatable :: ra
  ! variable for the structure
  integer :: nzone,isl,ill,nsl,nll
  integer, dimension(:), allocatable :: iphase
  integer :: ndc,vnp
  real(kind(0d0)) :: rmin,rmax
  real(kind(0d0)), dimension(:), allocatable :: vrmin,vrmax,qmu,qkappa,vra,rho,kappa
  real(kind(0d0)), dimension(:,:), allocatable :: rrho,vpv,vph,vsv,vsh,eta
  real(kind(0d0)), dimension(:), allocatable :: ecKx !3*Kx=3A-4N
  real(kind(0d0)), dimension(:), allocatable :: ecKy !3*Ky=3F+2N
  real(kind(0d0)), dimension(:), allocatable :: ecKz !3*Kz=2F+C
  real(kind(0d0)), dimension(:), allocatable :: mu,ecL,ecN,rhoinv,kappainv
  complex(kind(0d0)), dimension(:), allocatable :: coef1,coef2,coef
  ! variable for the periodic range
  integer :: np
  real(kind(0d0)) :: omega,omegai
  ! variable for the source
  integer :: spn,ns
  real(kind(0d0)):: spo,ecC0,ecF0,ecL0
  real(kind(0d0)), dimension(3,3) :: mt
  real(kind(0d0)), dimension(3) :: f
  complex(kind(0d0)), dimension(4) :: ya,yb,yc,yd
  ! variable for the matrix elements
  complex(kind(0d0)), dimension(:,:), allocatable :: a0,a1,a2,a,c
  real(kind(0d0)), dimension(:), allocatable :: t,h1,h2,h3,h4
  real(kind(0d0)), dimension(:), allocatable :: h1x,h1y,h1z,h2L,h2N,h3ax,h3ay,h3az,h4aL,h4aN
  real(kind(0d0)), dimension(:), allocatable :: h5ax,h5ay,h5az,h6aL,h6aN,h3x,h3y,h3z,h4L,h4N
  real(kind(0d0)), dimension(:), allocatable :: h5x,h5y,h5z,h6L,h6N,h7x,h7y,h7z,h8L,h8N
  real(kind(0d0)), dimension(:,:), allocatable :: h3mx,h3my,h3mz,h5mx,h5my,h5mz
  real(kind(0d0)), dimension(:,:), allocatable :: h4m1L,h4m1N,h4m2L,h4m2N,h6m1L,h6m1N,h6m2L,h6m2N
  real(kind(0d0)), dimension(:), allocatable :: p1,p2,p3
  complex(kind(0d0)), dimension(:), allocatable :: g0,d0
  complex(kind(0d0)), dimension(2) :: g0tmp,g0dertmp !forward
  ! variable for the stack point
  integer, dimension(:), allocatable :: isp,issp,ilsp,jssp,jsp,ksp,lsp
  integer :: isdr,jsdr,ildr,cista,cksta
  ! variables for the output stack point
  integer, dimension(:), allocatable :: istazone
  integer, dimension(:), allocatable :: ksta !output stack point for g
  integer, dimension(:), allocatable :: jsta !output stack point for d
  integer :: jjj
  ! variables for the gridding
  integer, dimension(:), allocatable :: jjdr,kkdr
  integer :: jdr,kdr
  real(kind(0d0)), dimension(:), allocatable :: vmin,gridpar,dzpar
  ! variables for l cut off
  integer :: kc,lsuf,sufzone,ismall
  real(kind(0d0)) :: maxamp
  ! variables for the numerical integration
  complex(kind(0d0)), dimension(4,4,10) :: anum,bnum
  ! other variables
  integer :: i,j,nn,lda,ier,itmp,jtmp,mtmp,kkdr0,nn0,ig2
  integer, dimension(12) :: ll,lli,llj
  real(kind(0d0)) :: eps,l2,lsq
  real(kind(0d0)), dimension(:), allocatable:: work,z,w,cwork
  !-----------------------------------------------------------------------
  complex(kind(0d0)), dimension(:,:,:), allocatable :: dvec,dvecdt,dvecdp,stress,displacement
  complex(kind(0e0)), dimension(:,:,:), allocatable :: stresssngl,displacementsngl
  real(kind(0d0)), dimension(:,:), allocatable :: plm
  real(kind(0d0)), dimension(:,:,:), allocatable :: plmtmp
  complex(kind(0d0)), dimension(3) :: u,udr,udt,udp
  complex(kind(0d0)), dimension(3,3) :: uder
  complex(kind(0d0)) :: rdvec(1:3,-2:2)
  data lda /4/
  data eps /-1.d0/

  include 'mpif.h'
  integer :: mysize,myrank,ierr
  integer :: worktag,deadtag,iwork,idone
  integer, dimension(MPI_STATUS_SIZE) :: mystat

  worktag = 1
  deadtag = 2

  call MPI_Init(ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD,mysize,ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD,myrank,ierr)

  if (myrank==0) then
    call date_and_time(datex,timex)
    write(*,'("start ",a4,"-",a2,"-",a2,"T",a2,":",a2,":",a4)') &
      datex(1:4),datex(5:6),datex(7:8),timex(1:2),timex(3:4),timex(5:8)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  write(*,'("proc ",i3," avail")') myrank
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myrank==0) then
    write(*,'("--> DSM : read input")')
    !!! read input
    call myinput(psvmodel,stationsinf,tlen,imin,imax,r0,r0lat,r0lon,i_source,i_geocentric)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_BCAST(psvmodel,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stationsinf,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(imin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(imax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(r0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(r0lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(r0lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(i_source,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(i_geocentric,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myrank==0) write(*,'("--> DSM : read model")')
  !!! read seismic model
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
  allocate(nlayer(nzone))
  allocate(iphase(nzone))
  allocate(vrmin(nzone))
  allocate(vrmax(nzone))
  allocate(rrho(4,nzone))
  allocate(vpv(4,nzone))
  allocate(vph(4,nzone))
  allocate(vsv(4,nzone))
  allocate(vsh(4,nzone))
  allocate(eta(4,nzone))
  allocate(qmu(nzone))
  allocate(qkappa(nzone))
  allocate(coef1(nzone))
  allocate(coef2(nzone))
  allocate(coef(nzone))
  allocate(jjdr(nzone))
  allocate(kkdr(nzone))
  allocate(vmin(nzone))
  allocate(gridpar(nzone))
  allocate(dzpar(nzone))
  allocate(isp(nzone))
  allocate(issp(nzone))
  allocate(ilsp(nzone))
  allocate(jssp(nzone))
  allocate(jsp(nzone)) 
  allocate(ksp(nzone))
  allocate(lsp(nzone))
  do i = 1, nzone
     read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), &
        vpv(1,i), vpv(2,i), vpv(3,i), vpv(4,i), vph(1,i), vph(2,i), vph(3,i), vph(4,i), &
        vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), &
        eta(1,i), eta(2,i), eta(3,i), eta(4,i), qmu(i), qkappa(i)
  enddo
  close(20)

  rmin = vrmin(1)
  rmax = vrmax(nzone)
  omegai = - dlog(1.d-2) / tlen
  
  if (myrank==0) write(*,'("--> DSM : identify model phase")')
  !!! identify model phase
  call calnl(nzone,vsv,iphase,nsl,nll)
  ndc = nzone - 1
  
  if (myrank==0) write(*,'("--> DSM : read stations")')
  !!! read stations file
  open (1,file=stationsinf,status='old',action='read',position='rewind')
  if(i_geocentric.eq.1) call translat (r0lat,r0lat)
  read(1,*)nsta
  r_n = nsta
  theta_n = nsta
  allocate(r_(r_n))
  allocate(A0sta(r_n))
  allocate(C0sta(r_n))
  allocate(F0sta(r_n))
  allocate(L0sta(r_n))
  allocate(N0sta(r_n))
  allocate(theta(theta_n))
  allocate(phi(theta_n))
  allocate(stla(nsta))
  allocate(stlo(nsta))
  allocate(updown(nsta))
  allocate(stress(6,6,nsta))
  allocate(displacement(3,6,nsta))
  allocate(stresssngl(6,6,nsta))
  allocate(displacementsngl(3,6,nsta))
  do i = 1,nsta
     read(1,*) r_(i),stla(i),stlo(i),updown(i)
     r_(i) = 6371.d0 -r_(i)
     if(i_geocentric.eq.1) call translat(stla(i),stla(i))
     !!! station coordinates relative to the source location
     call calthetaphi(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i))
     !!! ACFLN values at station
     call calstg4onedepth(nzone,vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,rmax,r_(i),updown(i), &
        A0sta(i),C0sta(i),F0sta(i),L0sta(i),N0sta(i))
  enddo
  close(1)

  !!! depths for stocking the Green function
  allocate(rrsta(3,r_n))
  allocate(iista(3,r_n))
  allocate(istazone(r_n))
  allocate(jsta(r_n))
  allocate(ksta(r_n))

  if (myrank==0) write(*,'("--> DSM : grid points")')
  !!! step in depth from omega and vmin
  call calgrid(nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,imax,1,tlen,vmin,gridpar,dzpar)
  !!! number of grid points by zone
  call calra_psv(nzone,vrmin,vrmax,iphase,dzpar,re,nnlayer,nslay,nllay,nlayer)
  allocate(ra(nnlayer+nzone+1))
  !!! calculate grid points location
  call calra2_psv(nnlayer,nzone,r_n,rmin,r0,nlayer,iphase,vrmin,vrmax,r_,cista,ciista,istazone,iista,ra,rrsta)
  nlay = nnlayer

  allocate(vra(nlay+2*nzone+1))
  allocate(rho(nlay+2*nzone+1))
  allocate(kappa(nlay+2*nzone+1))
  allocate(ecKx(nlay+2*nzone+1)) !3*Kx=3A-4N
  allocate(ecKy(nlay+2*nzone+1)) !3*Ky=3F+2N
  allocate(ecKz(nlay+2*nzone+1)) !3*Kz=2F+C
  allocate(mu(nlay+2*nzone+1))
  allocate(ecL(nlay+2*nzone+1))
  allocate(ecN(nlay+2*nzone+1))
  allocate(rhoinv(nlay+2*nzone+1))
  allocate(kappainv(nlay+2*nzone+1))
  allocate(a0(4,2*(2*(nslay+1)+(nllay+1)+2*nzone)))
  allocate(a1(4,2*(2*(nslay+1)+(nllay+1)+2*nzone)))
  allocate(a2(4,2*(2*(nslay+1)+(nllay+1)+2*nzone))) 
!!! GB GB
!  allocate(a(4,2*(nslay+1)+(nllay+1)))
!  allocate(c(2,(nslay+1)+(nllay+1)))
  allocate(a(4,2*(nslay+1)+(nllay+1)+nzone))
  allocate(c(2,(nslay+1)+(nllay+1)+nzone))
!!! GB GB
  allocate(t(8*nslay))
  allocate(h1x(8*nslay))
  allocate(h1y(8*nslay))
  allocate(h1z(8*nslay))
  allocate(h2L(8*nslay))
  allocate(h2N(8*nslay))
  allocate(h3ax(8*nslay))
  allocate(h3ay(8*nslay))
  allocate(h3az(8*nslay))
  allocate(h4aL(8*nslay))
  allocate(h4aN(8*nslay))
  allocate(h5ax(8*nslay))
  allocate(h5ay(8*nslay))
  allocate(h5az(8*nslay))
  allocate(h6aL(8*nslay))
  allocate(h6aN(8*nslay))
  allocate(h3x(8*nslay))
  allocate(h3y(8*nslay))
  allocate(h3z(8*nslay))
  allocate(h4L(8*nslay))
  allocate(h4N(8*nslay))
  allocate(h5x(8*nslay))
  allocate(h5y(8*nslay))
  allocate(h5z(8*nslay))
  allocate(h6L(8*nslay))
  allocate(h6N(8*nslay))
  allocate(h7x(8*nslay))
  allocate(h7y(8*nslay))
  allocate(h7z(8*nslay))
  allocate(h8L(8*nslay))
  allocate(h8N(8*nslay))
  allocate(h3mx(-2:1,2*(nslay+nzone)))
  allocate(h3my(-2:1,2*(nslay+nzone)))
  allocate(h3mz(-2:1,2*(nslay+nzone)))
  allocate(h5mx(-1:2,2*(nslay+nzone)))
  allocate(h5my(-1:2,2*(nslay+nzone)))
  allocate(h5mz(-1:2,2*(nslay+nzone)))
  allocate(h4m1L(-1:2,2*(nslay+nzone)))
  allocate(h4m1N(-1:2,2*(nslay+nzone)))
  allocate(h4m2L(-2:1,2*(nslay+nzone)))
  allocate(h4m2N(-2:1,2*(nslay+nzone)))
  allocate(h6m1L(-1:2,2*(nslay+nzone)))
  allocate(h6m1N(-1:2,2*(nslay+nzone)))
  allocate(h6m2L(-2:1,2*(nslay+nzone)))
  allocate(h6m2N(-2:1,2*(nslay+nzone)))
  allocate(p1(8*nllay))
  allocate(p2(8*nllay))
  allocate(p3(8*nllay))
  allocate(g0(2*(nslay+1)+(nllay+1)+nzone))
  allocate(d0((nslay+1)+(nllay+1)+nzone))
  allocate(work(8*nslay))
!!! GB GB
!  allocate(z(2*(nslay+1)+(nllay+1)))
!  allocate(w(2*(nslay+1)+(nllay+1)))
  allocate(z(2*(nslay+1)+(nllay+1)+nzone))
  allocate(w(2*(nslay+1)+(nllay+1)+nzone))
!!! GB GB
  allocate(cwork(4*(16*nslay+4*nllay)))

  !!! computing the stack points
  call calsp(nzone,ndc,nsl,nll,iphase,nlayer,nllay,isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr)
  !!! computing the source location
  call calspo(nlay,nzone,vrmax,iphase,nnlayer,ra,rmin,rmax,r0,isp,spo,spn)

  ! ******************* Computing the matrix elements *******************
  ! data initialization
  a = 0.d0
  t = 0.d0
  h1x = 0.d0
  h1y = 0.d0
  h1z = 0.d0
  h2L = 0.d0
  h2N = 0.d0
  h3ax = 0.d0
  h3ay = 0.d0
  h3az = 0.d0
  h4aL = 0.d0
  h4aN = 0.d0
  h5ax = 0.d0
  h5ay = 0.d0
  h5az = 0.d0
  h6aL = 0.d0
  h6aN = 0.d0
  h3x = 0.d0
  h3y = 0.d0
  h3z = 0.d0
  h4L = 0.d0
  h4N = 0.d0
  h5x = 0.d0
  h5y = 0.d0
  h5z = 0.d0
  h6L = 0.d0
  h6N = 0.d0
  h7x = 0.d0
  h7y = 0.d0
  h7z = 0.d0
  h8L = 0.d0
  h8N = 0.d0
  h3mx = 0.d0
  h3my = 0.d0
  h3mz = 0.d0
  h5mx = 0.d0
  h5my = 0.d0
  h5mz = 0.d0
  h4m1L = 0.d0
  h4m1N = 0.d0
  h4m2L = 0.d0
  h4m2N = 0.d0
  h6m1L = 0.d0
  h6m1N = 0.d0
  h6m2L = 0.d0
  h6m2N = 0.d0
  p1 = 0.d0
  p2 = 0.d0
  p3 = 0.d0

  ! computing the structure grid points
  call calstg(nlay,nzone,spn,rmax,r0,ra,nlayer,rrho,vpv,vph,vsv,vsh,eta,vra,rho,kappa,mu, &
    ecKx,ecKy,ecKz,ecL,ecN,ecC0,ecF0,ecL0,vnp)
  call calinv(nlay,nzone,vnp,rho,kappa,rhoinv,kappainv)

  isl = 0
  ill = 0
  do i = 1,ndc+1
     if (iphase(i)==1) then
        isl = isl+1
        itmp = isdr+issp(isl)
        call calmatc(nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)),t(itmp))
        call caltl(nlayer(i),vnp,vra,rho,ra(isp(i)),work(itmp))
        call calt(nlayer(i),t(itmp),work(itmp),t(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKx,0,0,0,ra(isp(i)),h1x(itmp))
        call calhl(nlayer(i),vnp,vra,ecKx,ra(isp(i)),work(itmp))
        call calt(nlayer(i),h1x(itmp),work(itmp),h1x(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKy,0,0,0,ra(isp(i)),h1y(itmp))
        call calhl(nlayer(i),vnp,vra,ecKy,ra(isp(i)),work(itmp))
        call calt(nlayer(i),h1y(itmp),work(itmp),h1y(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKz,0,0,0,ra(isp(i)),h1z(itmp))
        call calhl(nlayer(i),vnp,vra,ecKz,ra(isp(i)),work(itmp))
        call calt(nlayer(i),h1z(itmp),work(itmp),h1z(itmp))
        call calmatc(nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h2L(itmp))
        call calhl(nlayer(i),vnp,vra,ecL,ra(isp(i)),work(itmp))
        call calt(nlayer(i),h2L(itmp),work(itmp),h2L(itmp))
        call calmatc(nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h2N(itmp))
        call calhl(nlayer(i),vnp,vra,ecN,ra(isp(i)),work(itmp))
        call calt(nlayer(i),h2N(itmp),work(itmp),h2N(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKx,1,0,1,ra(isp(i)),h5ax(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKy,1,0,1,ra(isp(i)),h5ay(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKz,1,0,1,ra(isp(i)),h5az(itmp))
        call calmatc(nlayer(i),vnp,vra,ecL,1,0,1,ra(isp(i)),h6aL(itmp))
        call calmatc(nlayer(i),vnp,vra,ecN,1,0,1,ra(isp(i)),h6aN(itmp))
        call mtrnp(nlayer(i),h5ax(itmp),h3ax(itmp))
        call mtrnp(nlayer(i),h5ay(itmp),h3ay(itmp))
        call mtrnp(nlayer(i),h5az(itmp),h3az(itmp))
        call mtrnp(nlayer(i),h6aL(itmp),h4aL(itmp))
        call mtrnp(nlayer(i),h6aN(itmp),h4aN(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKx,2,1,1,ra(isp(i)),h7x(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKy,2,1,1,ra(isp(i)),h7y(itmp))
        call calmatc(nlayer(i),vnp,vra,ecKz,2,1,1,ra(isp(i)),h7z(itmp))
        call calmatc(nlayer(i),vnp,vra,ecL,2,1,1,ra(isp(i)),h8L(itmp))
        call calmatc(nlayer(i),vnp,vra,ecN,2,1,1,ra(isp(i)),h8N(itmp))
     else
        ill = ill+1
        itmp = ildr+ilsp(ill)
        call calmatc(nlayer(i),vnp,vra,rhoinv,2,1,1,ra(isp(i)),p1(itmp))
        call calmatc(nlayer(i),vnp,vra,rhoinv,0,0,0,ra(isp(i)),p2(itmp))
        call calhl(nlayer(i),vnp,vra,rhoinv,ra(isp(i)),work(itmp))
        call calt(nlayer(i),p2(itmp),work(itmp),p2(itmp))
        call calmatc(nlayer(i),vnp,vra,kappainv,2,0,0,ra(isp(i)),p3(itmp))
        call caltl(nlayer(i),vnp,vra,kappainv,ra(isp(i)),work(itmp))
        call calt(nlayer(i),p3(itmp),work(itmp),p3(itmp))
     endif
  enddo
  
  ! Computing the modified operator of the 1st derivative
  call caltstg(nlay,nzone,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,vra,kappa,ecKx,ecKy,ecKz,mu, &
    ecL,ecN)
  isl = 0
  do i = 1,ndc+1
     if (iphase(i)==1) then
        isl = isl+1
        itmp = isdr+issp(isl)
        jtmp = isp(i)+i-1
        call calh5(nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),work(itmp))
        call submat(nlayer(i),h5ax(itmp),work(itmp),h5x(itmp))
        call calh5(nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),work(itmp))
        call submat(nlayer(i),h5ay(itmp),work(itmp),h5y(itmp))
        call calh5(nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),work(itmp))
        call submat(nlayer(i),h5az(itmp),work(itmp),h5z(itmp))
        call calh5(nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),work(itmp))
        call submat(nlayer(i),h6aL(itmp),work(itmp),h6L(itmp))
        call calh5(nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),work(itmp))
        call submat(nlayer(i),h6aN(itmp),work(itmp),h6N(itmp))
        call mtrnp(nlayer(i),h5x(itmp),h3x(itmp))
        call mtrnp(nlayer(i),h5y(itmp),h3y(itmp))
        call mtrnp(nlayer(i),h5z(itmp),h3z(itmp))
        call mtrnp(nlayer(i),h6L(itmp),h4L(itmp))
        call mtrnp(nlayer(i),h6N(itmp),h4N(itmp))
        itmp = jsdr+jssp(isl)
        call calhm1(nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),h5mx(-1,itmp))
        call calhm1(nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),h5my(-1,itmp))
        call calhm1(nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),h5mz(-1,itmp))
        call calhm1(nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m1L(-1,itmp))
        call calhm1(nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m1N(-1,itmp))
        call calhm2(nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m2L(-2,itmp))
        call calhm2(nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m2N(-2,itmp))
        call mtrnp2(nlayer(i),1,2,h5mx(-1,itmp),h3mx(-2,itmp))
        call mtrnp2(nlayer(i),1,2,h5my(-1,itmp),h3my(-2,itmp))
        call mtrnp2(nlayer(i),1,2,h5mz(-1,itmp),h3mz(-2,itmp))
        call mtrnp2(nlayer(i),1,2,h6m1L(-1,itmp),h4m2L(-2,itmp))
        call mtrnp2(nlayer(i),1,2,h6m1N(-1,itmp),h4m2N(-2,itmp))
        call mtrnp2(nlayer(i),2,1,h6m2L(-2,itmp),h4m1L(-1,itmp))
        call mtrnp2(nlayer(i),2,1,h6m2N(-2,itmp),h4m1N(-1,itmp))
     endif
  enddo

allocate(plmtmp(1:3,0:3,1:theta_n))
allocate(plm(1:3,0:3))
allocate(dvec(1:3,-2:2,1:theta_n))
allocate(dvecdt(1:3,-2:2,1:theta_n))
allocate(dvecdp(1:3,-2:2,1:theta_n))

if (myrank==0) then
  do i = 1,mysize-1
    iwork = imax-(i-1) !imin+(i-1)
    call MPI_Send(iwork,1,MPI_INTEGER,i,worktag,MPI_COMM_WORLD,ierr)
    write(*,'("sent work ",i5," to proc ",i3," at ",e12.5)') iwork,i,dble(iwork)/tlen
  enddo
  do i = mysize,(imax-imin)+1
    iwork = imax-(i-1) !imin+(i-1)
    call MPI_Recv(idone,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mystat,ierr)
    call MPI_Send(iwork,1,MPI_INTEGER,mystat(MPI_SOURCE),worktag,MPI_COMM_WORLD,ierr)
    write(*,'("sent work ",i5," to proc ",i3," at ",e12.5)') iwork,mystat(MPI_SOURCE),dble(iwork)/tlen
  enddo
  do i = 1,mysize-1
    call MPI_Recv(idone,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mystat,ierr)
    call MPI_Send(0,0,MPI_INTEGER,mystat(MPI_SOURCE),deadtag,MPI_COMM_WORLD,ierr)
    write(*,'("stop proc ",i3)') mystat(MPI_SOURCE)
  enddo
else

  coutfile = "out_stress_PSV.bin"
  open(1,file=coutfile,status='unknown',form='unformatted',access='direct',recl=6*6*nsta*2*kind(0e0),iostat=ierr)
  coutfile = "out_displ_PSV.bin"
  open(2,file=coutfile,status='unknown',form='unformatted',access='direct',recl=3*6*nsta*2*kind(0e0),iostat=ierr)

  do !!! worker loop
    call MPI_Recv(iwork,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,mystat,ierr)
    if (mystat(MPI_TAG) == deadtag) exit
    write(*,'("proc ",i3," working on ",i5," at ",e12.5)') myrank,iwork,dble(i)/tlen
    i = iwork

    stress = cmplx(0.d0)
    displacement = cmplx(0.d0)
    stresssngl = cmplx(0.e0)
    displacementsngl = cmplx(0.e0)
     
    omega = 2.d0*pi*dble(i)/tlen

    if (i/=0) then       
      call callsuf(omega,nzone,vrmax,vsv,lsuf)
      call calcoef(nzone,omega,qmu,qkappa,coef1,coef2,coef)  
      mtmp = isp(spn)+int(spo)
      if (spo==int(spo)) mtmp = mtmp-1
      call calabnum(omega,omegai,rmax,rrho(1,spn),vpv(1,spn),vph(1,spn),vsv(1,spn),vsh(1,spn), &
        eta(1,spn),ra(mtmp),r0,coef1(spn),coef2(spn),anum(1,1,1),bnum(1,1,1))
      !!! computing the matrix elements independent of l
      isl = 0
      ill = 0
      do j = 1,ndc+1
        if (iphase(j)==1) then
          isl = isl+1
          itmp = isdr+issp(isl)
          jtmp = jdr+jsp(j)
          mtmp = kdr+ksp(j)
          call cala0(nlayer(j),omega,omegai,t(itmp),h1x(itmp),h1y(itmp),h1z(itmp),h2L(itmp), &
            h2N(itmp),h3ax(itmp),h3ay(itmp),h3az(itmp),h4aL(itmp),h4aN(itmp),h5ax(itmp), &
            h5ay(itmp),h5az(itmp),h6aL(itmp),h6aN(itmp),h7x(itmp),h7y(itmp),h7z(itmp),h8L(itmp), &
            h8N(itmp),coef1(j),coef2(j),cwork(jtmp))
          call overlapa(nlayer(j),cwork(jtmp),a0(1,mtmp))
          call cala1(nlayer(j),h1x(itmp),h1y(itmp),h1z(itmp),h2L(itmp),h2N(itmp),h3x(itmp), &
            h3y(itmp),h3z(itmp),h4L(itmp),h4N(itmp),h5x(itmp),h5y(itmp),h5z(itmp),h6L(itmp), &
            h6N(itmp),coef1(j),coef2(j),cwork(jtmp))
          call overlapa(nlayer(j),cwork(jtmp),a1(1,mtmp))
          call cala2(nlayer(j),h1x(itmp),h1y(itmp),h1z(itmp),h2L(itmp),h2N(itmp),coef1(j), &
            coef2(j),cwork(jtmp))
          call overlapa(nlayer(j),cwork(jtmp),a2(1,mtmp))
              jtmp = jsdr+jssp(isl)
          call calhml(nlayer(j),coef1(j),coef2(j),h3mx(-2,jtmp),h3my(-2,jtmp),h3mz(-2,jtmp), &
            h5mx(-1,jtmp),h5my(-1,jtmp),h5mz(-1,jtmp),h4m1L(-1,jtmp),h4m1N(-1,jtmp), &
            h4m2L(-2,jtmp),h4m2N(-2,jtmp),h6m1L(-1,jtmp),h6m1N(-1,jtmp),h6m2L(-2,jtmp), &
            h6m2N(-2,jtmp),a1(1,mtmp))
        else
          ill = ill+1
          itmp = ildr+ilsp(ill)
          jtmp = jdr+jsp(j)
          mtmp = kdr+ksp(j)
          call calb0(nlayer(j),omega,omegai,p1(itmp),p3(itmp),coef(j),cwork(jtmp))
          call overlapb(nlayer(j),cwork(jtmp),a0(1,mtmp))
          call calb2(nlayer(j),omega,omegai,p2(itmp),coef(j),cwork(jtmp))
          call overlapb(nlayer(j),cwork(jtmp),a2(1,mtmp))
        endif
      enddo
    
      kc = 1
      ismall = 0
      maxamp = -1.d0
      plm = 0.
      plmtmp = 0.
      do l=0,maxlmax    ! l-loop start

        do itheta = 1,theta_n
          call clPLM_modified(l,(theta(itheta)/180.d0*pi),plmtmp(1:3,0:3,itheta),plm(1:3,0:3))
          call caldvec_dejaplm(l,(theta(itheta)/180.d0*pi),(phi(itheta)/180.d0*pi),plm(1:3,0:3), &
                dvec(1:3,-2:2,itheta),dvecdt(1:3,-2:2,itheta),dvecdp(1:3,-2:2,itheta))
        enddo

        l2 = dble(l)*dble(l+1)
        lsq = dsqrt( l2 )

        rdvec = cmplx(0.d0)
        call caldveczero(l,rdvec(1:3,-2:2))

        ! computing the coefficient matrix elements
        ! --- renewing  mdr
        if ( mod(l,50).eq.0 )  then
          call calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
          call calspdr(nzone,nzone,iphase,nlayer,jjdr,kkdr )
          do ir_=1,r_n
            ksta(ir_) = kkdr(istazone(ir_))+2*iista(1,ir_) - 1
          enddo
          cksta = kkdr(istazone(cista))+2*iista(1,cista) - 1
          nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
        endif

        !     computing the matrix elements
        call cala( nzone,ndc,iphase,nlayer,kkdr,kdr,ksp,l2,lsq,nn,a0,a1,a2,a )
        ! computing the boundary condition elements
        call calbc( nzone,ndc,vrmax,iphase,kkdr,a )
           
        jtmp = kkdr(spn) + 2 * int(spo)
        mtmp = isp(spn) + int(spo)
        if ( spo.eq.int(spo) ) then
          jtmp = jtmp - 2
          mtmp = mtmp - 1
        endif
  
        call calya( anum(1,1,1),bnum(1,1,1),l2,ra(mtmp),r0,ya,yb,yc,yd )
        do m=-2,2        ! m-loop start
          if ( iabs(m).le.iabs(l) ) then
            ig2 = 0
            if ( l.eq.0 ) then ! l-branch for calu (l=0) 
              !  rearranging the matrix elements 
              do imt = 1,6
                g0 = cmplx(0.d0)
                itmp = 1
                if (rmin==0.d0) itmp = 2

                if (i_source==0) then
                  !!! moment tensor source
                  call setmt(imt,mt)
                  call calg(l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0,mt,g0(jtmp))
                  !if (iwork==100.and.l<10) write(*,'(3i6,10e12.4)') l,m,imt,g0(jtmp),g0(jtmp+1),g0(jtmp+2),g0(jtmp+3)
                  call rea2(nn,a,g0,c,d0,nzone,iphase,kkdr,spn,kkdr0,nn0,r_n,r_n,istazone,iista,jsta)
                  call dcsymbdl0(c(1,itmp),1,nn0-itmp+1,1,eps,z(itmp),w(itmp),ll,lli,llj,ier)
                  if ((abs(m)==0.and.(imt==1.or.imt==2.or.imt==3)).or. &
                    (abs(m)==1.and.(imt==4.or.imt==5)).or. &
                    (abs(m)==2.and.(imt==2.or.imt==3.or.imt==6))) then
                    call mydcsbdlv1(c(1,itmp),d0(itmp),nn0-itmp+1,z(itmp))
                  endif
                elseif (i_source==1) then
                  !!! force vector source
                  if (imt>3.or.iabs(m)==2) cycle
                  call setf(imt,f)
                  call calg_force(l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0,f,g0(jtmp))
                  !if (iwork==100.and.l<10) write(*,'(3i6,10e12.4)') l,m,imt,g0(jtmp),g0(jtmp+1),g0(jtmp+2),g0(jtmp+3)
                  call rea2(nn,a,g0,c,d0,nzone,iphase,kkdr,spn,kkdr0,nn0,r_n,r_n,istazone,iista,jsta)
                  call dcsymbdl0(c(1,itmp),1,nn0-itmp+1,1,eps,z(itmp),w(itmp),ll,lli,llj,ier)
                  if ((abs(m)==0.and.imt==1).or.(abs(m)==1.and.(imt==2.or.imt==3))) then
                    call mydcsbdlv1(c(1,itmp),d0(itmp),nn0-itmp+1,z(itmp))
                  endif
                endif

                do ir_=1,r_n
                  g0tmp = cmplx(0.d0)
                  g0dertmp = cmplx(0.d0)
                  call interpolate( 1,0,r_(ir_),rrsta(1,ir_),d0(jsta(ir_)),g0tmp(1))
                  call interpolate( 1,1,r_(ir_),rrsta(1,ir_),d0(jsta(ir_)),g0dertmp(1))
                  itheta = ir_
                  u = cmplx(0.d0)
                  udr = cmplx(0.d0)
                  udt = cmplx(0.d0)
                  udp = cmplx(0.d0)
                  uder = cmplx(0.d0)
                  call calup0(g0tmp(1),dvec(1:3,m,itheta),u(1:3))
                  call calup0(g0dertmp(1),dvec(1:3,m,itheta),udr(1:3))
                  call calup0(g0tmp(1),dvecdt(1:3,m,itheta),udt(1:3))
                  call calup0(g0tmp(1),dvecdp(1:3,m,itheta),udp(1:3))
                  call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3), &
                    r_(ir_),theta(itheta)/180.d0*pi)
                  call udertoStress(uder(1:3,1:3),stress(1:6,imt,itheta),A0sta(itheta),C0sta(itheta), &
                    F0sta(itheta),L0sta(itheta),N0sta(itheta))
                  displacement(1:3,imt,itheta) = u(1:3)+displacement(1:3,imt,itheta)
                enddo
              enddo ! imt-loop
            else ! for l!=0
              do imt = 1,6
                g0 = cmplx(0.d0)
                itmp = 1
                if (rmin==0.d0) itmp = 3

                if (i_source==0) then
                  !!! moment tensor source
                  call setmt(imt,mt)
                  call calg(l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0,mt,g0(jtmp))
                  ! computing forward propagating component (l!=0)                              
                  if ((m==-2.or.m==-l).and.(ig2==0)) then
                    call dcsymbdl0(a(1,itmp),3,nn-itmp+1,6,eps,z(itmp),w(itmp),ll,lli,llj,ier)
                    ig2 = 1
                  endif
                  if ((abs(m)==0.and.(imt==1.or.imt==2.or.imt==3)).or. &
                      (abs(m)==1.and.(imt==4.or.imt==5)).or. &
                      (abs(m)==2.and.(imt==2.or.imt==3.or.imt==6))) then
                    call mydcsbdlv3(a(1,itmp),g0(itmp),nn-itmp+1,z(itmp))
                  endif
                  if (imt==1) call calamp(g0(ksta(r_n)-1),l,lsuf,maxamp,ismall,ratl)
                elseif (i_source==1) then
                  !!! force vector source
                  if (imt>3.or.iabs(m)==2) cycle
                  call setf(imt,f)
                  call calg_force(l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0,f,g0(jtmp))
                  if (m==-1.and.ig2==0) then
                    call dcsymbdl0(a(1,itmp),3,nn-itmp+1,6,eps,z(itmp),w(itmp),ll,lli,llj,ier)
                    ig2 = 1
                  endif
                  if ((abs(m)==0.and.imt==1).or.(abs(m)==1.and.(imt==2.or.imt==3))) then
                    call mydcsbdlv3(a(1,itmp),g0(itmp),nn-itmp+1,z(itmp))
                  endif
                endif

                do ir_=1,r_n ! stack point
                  g0tmp = cmplx(0.d0)
                  g0dertmp = cmplx(0.d0)
                  call interpolate(2,0,r_(ir_),rrsta(1,ir_),g0(ksta(ir_)-1),g0tmp(1:2))
                  call interpolate(2,1,r_(ir_),rrsta(1,ir_),g0(ksta(ir_)-1),g0dertmp(1:2) )
                  itheta = ir_
                  u = cmplx(0.d0)
                  udr = cmplx(0.d0)
                  udt = cmplx(0.d0)
                  udp = cmplx(0.d0)
                  uder = cmplx(0.d0)
                  call calup(g0tmp(1),g0tmp(2),lsq,dvec(1:3,m,itheta),u(1:3))
                  call calup(g0dertmp(1),g0dertmp(2),lsq,dvec(1:3,m,itheta),udr(1:3))
                  call calup(g0tmp(1),g0tmp(2),lsq,dvecdt(1:3,m,itheta),udt(1:3))
                  call calup(g0tmp(1),g0tmp(2),lsq,dvecdp(1:3,m,itheta),udp(1:3))
                  call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3), &
                    r_(ir_),theta(itheta)/180.d0*pi)
                  call udertoStress(uder(1:3,1:3),stress(1:6,imt,itheta),A0sta(itheta),C0sta(itheta), &
                    F0sta(itheta),L0sta(itheta), N0sta(itheta))
                  displacement(1:3,imt,itheta) = u(1:3)+displacement(1:3,imt,itheta)
                enddo ! stack point
              enddo ! mt-loop
            endif ! l-branch for calu
          endif
        enddo ! m-loop end
      enddo ! l-loop end        


      stress(1:6,1:6,1:nsta) = stress(1:6,1:6,1:nsta)/cmplx(0,omega) 
      stresssngl(1:6,1:6,1:nsta) = stress(1:6,1:6,1:nsta)
      displacementsngl(1:3,1:6,1:nsta) = displacement(1:3,1:6,1:nsta)

    endif ! i/=0

    write(1,rec=i+1,iostat=ierr) stresssngl(1:6,1:6,1:nsta)
    write(2,rec=i+1,iostat=ierr) displacementsngl(1:3,1:6,1:nsta)

    idone = iwork
    call MPI_Send(idone,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

  enddo ! frequencies loop

  close(1,iostat=ierr)
  close(2,iostat=ierr)

endif ! worker/slave

if (myrank==0) then
  call date_and_time(datex,timex)
  write(*,'("end ",a4,"-",a2,"-",a2,"T",a2,":",a2,":",a4)') datex(1:4),datex(5:6),datex(7:8), &
        timex(1:2),timex(3:4),timex(5:8)
endif

call MPI_Finalize(ierr)

stop
 
end program myTraPSV
