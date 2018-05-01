program myTraSH


!-----------------------------------------------------------------------
!     
!
!  
!
!       
!                                               2002.10.KAWAI Kenji
!                                               2009.6. FUJI Nobuaki
!                                               2011.9. FUJI Nobuaki
!                                               
!
!
!                 
!
!-----------------------------------------------------------------------

  implicit none

!-------------------------<< input matrix >>----------------------------------

  character(120) :: psvmodel,stationsinf
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)), parameter :: re=1.d-2, ratc=1.d-10, ratl=1.d-4
  integer, parameter :: maxlmax = 20000
 
  real(kind(0d0)) :: tlen 
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: r0lat, r0lon
  real(kind(0d0)), allocatable :: stla(:),stlo(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:)
  real(kind(0d0)), allocatable :: A0sta(:),C0sta(:),F0sta(:),L0sta(:),N0sta(:)
  integer, allocatable :: updown(:) 
  real(kind(0d0)), allocatable :: rrsta(:,:)
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n,ciista, ir_,ir0,imt,icomp,idepth,itheta,theta_n,nsta
  
  character(120) :: coutfile
  integer :: imin,imax
  integer :: i_source,i_geocentric
  integer :: kpa
  integer :: nzone
  integer :: iimax,ii

  integer :: i ,j, ier,jj
  real(kind(0d0)) :: dummy
  real(kind(0d0)), allocatable :: vrmin(:), vrmax(:)
  real(kind(0d0)), allocatable :: rrho(:,:), vsv(:,:), vsh(:,:), qmu(:)
  real(kind(0d0)), allocatable :: vra (:), rho(:), ecL(:), ecN(:)
  real(kind(0d0)), allocatable :: gvra(:,:), grho(:,:), gecL(:,:), gecN(:,:),gra(:,:)
  complex(kind(0d0)), allocatable :: coef(:),cwork(:)
  real(kind(0d0)) :: rmin, rmax
  real(kind(0d0)), allocatable :: vmin(:), gridpar(:), dzpar(:)
  complex(kind(0d0)), allocatable :: tmpc(:)
  real(kind(0d0)) :: maxamp
  real(kind(0d0)) :: omegai
  real(kind(0d0)), allocatable :: ra(:)
  integer :: nnlayer, vnp,nn
  integer, allocatable :: nlayer(:), iphase(:)
  integer :: ioutercore

  ! variables pour des points stackes 
  integer, allocatable :: isp(:),jsp(:),ins(:)

  ! variables pour la source
 
  integer, allocatable :: spn(:),ns(:)
  real(kind(0d0)) :: mt(3,3),f(3),lsq
  real(kind(0d0)), allocatable :: mu0(:),spo(:)

!-----------------------------------------------------------------------
  ! variables pour des elements de matrice 
  complex(kind(0d0)), allocatable :: a0(:,:), a2(:,:), a(:,:),dr(:),z(:)
  real(kind(0d0)), allocatable :: t(:), h1(:), h2(:), h3(:), h4(:), work(:)
  real(kind(0d0)), allocatable :: gt(:,:),gh1(:,:),gh2(:,:),gh3(:,:),gh4(:,:)
  complex(kind(0d0)),allocatable :: aa(:,:), ga(:,:),ga2(:,:,:),gdr(:,:)
  complex(kind(0d0)), allocatable :: g0(:)
  complex(kind(0d0)) :: g0tmp, g0dertmp
  ! la frequence
  real(kind(0d0)) :: omega
  integer :: lsuf

  ! des autres 
  integer :: lda 
  integer :: kc, ismall,m, l,ig2
  real(kind(0d0)) :: eps 
  
!-----------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: bvec(:,:,:),bvecdt(:,:,:),bvecdp(:,:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:)
  complex(kind(0d0)), allocatable :: stress(:,:,:), displacement(:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl(:,:,:), displacementsngl(:,:,:)
  complex(kind(0d0))::u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3),rvec(1:3,-2:2)

  data lda / 2 /
  data eps / -1.d0 /

  integer :: iwork

    call date_and_time(datex,timex)
    write(*,'("start ",a4,"-",a2,"-",a2,"T",a2,":",a2,":",a4)') datex(1:4),datex(5:6),datex(7:8), &
        timex(1:2),timex(3:4),timex(5:8)

    write(*,'("--> DSM : read input")')
    call myinput(psvmodel,stationsinf,tlen,imin,imax,r0min,r0lat,r0lon,i_source,i_geocentric)
    r0max = r0min
    r0delta = 20.d0

  write(*,'("--> DSM : read model")')
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))
  do i = 1, nzone
     read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), &
        dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
        vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), &
        dummy, dummy, dummy, dummy, qmu(i), dummy
     if ((vsv(1,i).eq.0.d0).and.(vsv(2,i).eq.0.d0).and.(vsv(3,i).eq.0.d0).and.(vsv(4,i).eq.0.d0)) then
        iphase(i) = 2
        ioutercore = i
     else
        iphase(i) = 1
     endif
  enddo
  close(20)

  ! CAUTION: this program can only calculate for solid media (SH) for this moment
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
  nzone = nzone - ioutercore
  deallocate(vrmin,vrmax,rrho,vsv,vsh,qmu,vmin,gridpar,dzpar,nlayer,iphase,isp,jsp,coef)
  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))
  do i = 1,ioutercore
     read (20, *) dummy, dummy, dummy, dummy, dummy, dummy, &
        dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
        dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
        dummy, dummy, dummy, dummy, dummy, dummy
  enddo
  do i = 1, nzone
     read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), &
        dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
        vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), &
        dummy, dummy, dummy, dummy, qmu(i), dummy
  enddo
  close(20)

  rmin = vrmin(1)
  rmax = vrmax(nzone)
  omegai = - dlog(1.d-2) / tlen

  write(*,'("--> DSM : read stations")')
  open (1,file=stationsinf,status='old',action='read',position='rewind')
  if(i_geocentric.eq.1) call translat (r0lat,r0lat)
  read(1,*)nsta
  r_n = nsta
  theta_n = nsta
  allocate(r_(1:r_n))
  allocate(A0sta(1:r_n))
  allocate(C0sta(1:r_n))
  allocate(F0sta(1:r_n))
  allocate(L0sta(1:r_n))
  allocate(N0sta(1:r_n))
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stla(1:nsta))
  allocate(stlo(1:nsta))
  allocate(updown(1:nsta))
  allocate(stress(1:6,1:6,1:nsta))
  allocate(displacement(1:3,1:6,1:nsta))
  allocate(stresssngl(1:6,1:6,1:nsta))
  allocate(displacementsngl(1:3,1:6,1:nsta))
  do i = 1,nsta
     read(1,*) r_(i),stla(i),stlo(i),updown(i)
     r_(i) = 6371.d0 -r_(i)
     if(i_geocentric.eq.1) call translat(stla(i),stla(i))
     call calthetaphi(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i))
     call calstg4onedepth(nzone,nzone,nzone,vrmin,vrmax,iphase,rrho,0.d0,0.d0,vsv,vsh,1.d0,rmax,r_(i),updown(i), &
        A0sta(i),C0sta(i),F0sta(i),L0sta(i),N0sta(i))
  enddo
  close(1)
  
  allocate(rrsta(1:3,1:r_n))
  allocate(iista(1:3,1:r_n))
  allocate(bvec(1:3,-2:2,1:theta_n))
  allocate(bvecdt(1:3,-2:2,1:theta_n))
  allocate(bvecdp(1:3,-2:2,1:theta_n))
  allocate(plm(1:3,0:3,1:theta_n))

  ! source depths

  r0_n =  int((r0max-r0min)/r0delta)+1
  
  allocate(r0(1:r0_n))
  allocate(spo(1:r0_n))
  allocate(spn(1:r0_n))
  allocate(ns(1:r0_n))
  allocate(mu0(1:r0_n))
  allocate(ins(1:r0_n))
  allocate(gra(1:3,1:r0_n)) 
  allocate(gvra(1:3,1:r0_n)) 
  allocate(grho(1:3,1:r0_n))
  allocate(gecL(1:3,1:r0_n))
  allocate(gecN(1:3,1:r0_n)) 
  allocate(gt(1:8,1:r0_n))
  allocate(gh1(1:8,1:r0_n))
  allocate(gh2(1:8,1:r0_n))
  allocate(gh3(1:8,1:r0_n))
  allocate(gh4(1:8,1:r0_n)) 
  allocate(aa(1:4,1:r0_n))
  allocate(ga(1:8,1:r0_n))
  allocate(ga2(1:2,1:3,1:r0_n))
  allocate(gdr(1:3,r0_n))

  
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo

  ! computation de nombre et la location des points de grid

  write(*,'("--> DSM : compute grid points")')

  iimax = imax

  call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,iimax,1,tlen,vmin,gridpar,dzpar )
  call calra(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,re )
  allocate(ra(1:nnlayer+nzone+1))
  call calra2(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,ra,re,r_n,r_,rrsta, iista)
  ! computation de points stackes et la location de la source
  call calsp( nzone-1,nlayer,isp,jsp )
  do ir0 = 1,r0_n
     call calspo( nzone-1,vrmax,nnlayer,r0(ir0),rmin,rmax,ra,isp,spo(ir0),spn(ir0) )
     call calgra( isp,ra,r0(ir0),spn(ir0),spo(ir0),gra(1:3,ir0))
  enddo

  ! computation des elements de matrice

  write(*,'("--> DSM : compute matrix elements")')
  allocate(vra(1:nnlayer+2*nzone+1))
  allocate(rho(1:nnlayer+2*nzone+1))
  allocate(ecL(1:nnlayer+2*nzone+1))
  allocate(ecN(1:nnlayer+2*nzone+1))
  allocate(a0(1:2,1:nnlayer+1))
  allocate(a2(1:2,1:nnlayer+1))
  allocate(a(1:2,1:nnlayer+1))
  allocate(t(1:4*nnlayer))
  allocate(cwork(1:4*nnlayer))
  allocate(h1(1:4*nnlayer))
  allocate(h2(1:4*nnlayer))
  allocate(h3(1:4*nnlayer))
  allocate(h4(1:4*nnlayer))
  allocate(work(1:4*nnlayer))
  allocate(tmpc(nnlayer+1))
  allocate(g0(1:nnlayer+1))
  allocate(dr(1:nnlayer+1))
  allocate(z(1:nnlayer+1))  
  call calstg( nzone,rrho,vsv,vsh,nnlayer,nlayer,ra,rmax,vnp,vra,rho,ecL,ecN)
  do ir0 = 1, r0_n
     call calgstg(nzone,nnlayer,spn(ir0),rrho,vsv,vsh,gra(1:3,ir0),gvra(1:3,ir0),rmax, &
        grho(1:3,ir0),gecL(1:3,ir0),gecN(1:3,ir0),r0(ir0),mu0(ir0)) 
  enddo

  do i= 1, nzone
     call calmatc(nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)),t(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,2,1,1,ra(isp(i)),h1(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,1,1,0,ra(isp(i)),h2(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h3(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h4(jsp(i)),work(jsp(i)))
     call caltl(nlayer(i),vnp,vra,rho,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),t(jsp(i)),work(jsp(i)),t(jsp(i)))
     call calhl(nlayer(i),vnp,vra,ecL,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),h3(jsp(i)),work(jsp(i)),h3(jsp(i)))
     call calhl(nlayer(i),vnp,vra,ecN,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),h4(jsp(i)),work(jsp(i)),h4(jsp(i)))
  enddo
  do ir0 = 1, r0_n
     call calmatc( 2,3,gvra(1:3,ir0),grho(1:3,ir0),2,0,0,gra(1:3,ir0),gt(1:8,ir0), work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),2,1,1,gra(1:3,ir0),gh1(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),1,1,0,gra(1:3,ir0),gh2(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),0,0,0,gra(1:3,ir0),gh3(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecN(1:3,ir0),0,0,0,gra(1:3,ir0),gh4(1:8,ir0),work )
     call caltl(2,3,gvra(1:3,ir0),grho(1:3,ir0),gra(1:3,ir0),work )
     call calt(2,gt(1:8,ir0),work,gt(1:8,ir0))
     call calhl(2,3,gvra(1:3,ir0),gecL(1:3,ir0),gra(1:3,ir0),work)
     call calt( 2,gh3(1:8,ir0),work,gh3(1:8,ir0))
     call calhl(2,3,gvra(1:3,ir0),gecN(1:3,ir0),gra(1:3,ir0),work )
     call calt( 2,gh4(1:8,ir0),work,gh4(1:8,ir0))  
  enddo


  !computation de la dislocation
  nn = nnlayer + 1
  do ir0 = 1, r0_n
     ns(ir0) = isp(spn(ir0)) + dint(spo(ir0))
     ins(ir0) = 4 * ns(ir0) - 3
  enddo

  do iwork = imax,imin,-1 !!! worker loop
    write(*,'("proc ",i3," working on ",i5," at ",e12.5)') 0,iwork,dble(iwork)/tlen
    i = iwork

     stress = cmplx(0.d0)
     displacement = cmplx(0.d0)
     stresssngl=cmplx(0.e0)
     displacementsngl=cmplx(0.e0)
     
     omega = 2.d0 * pi * dble(i)/tlen
     if ( i.ne.0 ) then
        call callsuf(omega,nzone,vrmax,vsv,lsuf)       
        call calcoef( nzone,omega,qmu,coef)
        plm = 0.d0
        a0 = 0.d0
        a2 = 0.d0
        do j = 1, nzone
           call cala0( nlayer(j),omega,omegai,t(jsp(j)), h1(jsp(j)),&
                & h2(jsp(j)), h3(jsp(j)), h4(jsp(j)),coef(j), cwork(jsp(j)) )
           call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
           call cala2( nlayer(j),h4(jsp(j)),coef(j), cwork(jsp(j)) )
           call overlap( nlayer(j),cwork(jsp(j)),a2( 1,isp(j) ) )
        enddo

        kc = 1
        ismall = 0
        maxamp = -1.d0
        do l = 0, maxlmax ! l-loop commence
           lsq = dsqrt(dble(l)*dble(l+1))
           do itheta = 1,theta_n
              call calbvec(l,(theta(itheta)/180.d0*pi),(phi(itheta)/180.d0*pi),plm(1,0,itheta),bvec(1,-2,itheta), &
                bvecdt(1,-2,itheta),bvecdp(1,-2,itheta))
           enddo

           rvec = cmplx(0.d0)
           !!! GB GB call calbveczero(l,rvec(1,-2))
           call calbveczero(l,rvec(1,-2))
           !if (i==100) write(*,'(i8,10e10.2)') l,rvec(3,-2:2)
           !!!! GB GB
           tmpc = 0.d0
           a = 0.d0
           ga2 = 0.d0
           
           call cala( nn,l,lda,a0,a2,a )
           do ir0 = 1, r0_n
              call calga( 1,omega,omegai,l,t(ins(ir0)),h1(ins(ir0)),h2(ins(ir0)),h3(ins(ir0)),h4(ins(ir0)), &
                coef(spn(ir0)),aa(1:4,ir0))
              call calga( 2,omega,omegai,l,gt(1:8,ir0),gh1(1:8,ir0),gh2(1:8,ir0),gh3(1:8,ir0),gh4(1:8,ir0), &
                coef(spn(ir0)),ga(1:8,ir0))
              call overlap( 2,ga(1:8,ir0),ga2(1:2,1:3,ir0))
           enddo
           
           do m = -2, 2 ! m-loop commence       
              if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
                 do ir0 = 1,r0_n
                    ig2 = 0
                    do imt = 2,6
                       g0 = cmplx(0.d0)
                       if (i_source==0) then
                         !!! moment tensor source
                         call setmt(imt,mt)
                         call calg2(l,m,spo(ir0),r0(ir0),mt,mu0(ir0),coef(spn(ir0)),ga(1:8,ir0),aa(1:4,ir0), &
                           ga2(1:2,1:3,ir0),gdr(1:3,ir0),g0(isp(spn(ir0))),ig2,i)
                         !if (iwork==100.and.l<10) write(*,'(3i6,4e12.4)') l,m,imt,g0(isp(spn(ir0))),g0(isp(spn(ir0))+1)
                         if ((m==-2).or.(m==-l)) then
                           if (ig2==1) then
                             call dclisb0(a,nn,1,lda,g0,eps,dr,z,ier)
                             ig2 = ig2+1
                           else
                             call dcsbsub0(a,nn,1,lda,g0,eps,dr,z,ier)
                           endif
                         else
                           call dcsbsub0(a,nn,1,lda,g0,eps,dr,z,ier)
                         endif
                         if ((imt==3).and.(ir0==r0_n)) call calamp(g0(nn-1),l,lsuf,maxamp,ismall,ratl)
                       elseif (i_source==1) then
                         !!! force vector source
                         if (imt>3.or.iabs(m)==2) cycle
                         call setf(imt,f)
                         call calg_force(l,m,spo(ir0),r0(ir0),f,mu0(ir0),coef(spn(ir0)),ga(1:8,ir0),aa(1:4,ir0), &
                           ga2(1:2,1:3,ir0),gdr(1:3,ir0),g0(isp(spn(ir0))),ig2,i)
                         !if (iwork==100.and.l<10) write(*,'(3i6,4e12.4)') l,m,imt,g0(isp(spn(ir0))),g0(isp(spn(ir0))+1)
                         if (m==-1.and.ig2==1) then
                           call dclisb0(a,nn,1,lda,g0,eps,dr,z,ier)
                           ig2 = ig2+1
                         else
                           call dcsbsub0(a,nn,1,lda,g0,eps,dr,z,ier)
                         endif
                       endif

                       do ir_= 1,r_n                    
                          g0tmp = 0.d0
                          g0dertmp = 0.d0
                          call interpolate(1,0,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0tmp)
                          call interpolate(1,1,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0dertmp)
                          itheta = ir_
                          u = cmplx(0.d0)
                          udr = cmplx(0.d0)
                          udt = cmplx(0.d0)
                          udp = cmplx(0.d0)
                          uder = cmplx(0.d0)
                          call calu(g0tmp,lsq,bvec(1:3,m,itheta),u(1:3))
                          call calu(g0dertmp,lsq,bvec(1:3,m,itheta),udr(1:3))
                          call calu(g0tmp,lsq,bvecdt(1:3,m,itheta),udt(1:3))
                          call calu(g0tmp,lsq,bvecdp(1:3,m,itheta),udp(1:3))
                          call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3), &
                            r_(ir_),theta(itheta)/180.d0*pi)
                          call udertoStress(uder(1:3,1:3),stress(1:6,imt,itheta),A0sta(itheta),C0sta(itheta), &
                            F0sta(itheta),L0sta(itheta),N0sta(itheta))
                          displacement(1:3,imt,itheta) = u(1:3)+displacement(1:3,imt,itheta)
                       enddo ! ir_-loop termine
                    enddo ! imt-loop termine
                 enddo !ir0-loop termine     
              endif
           enddo ! m-loop termine              
        enddo !l-loop termine                       

        stress(1:6,1:6,1:nsta) = stress(1:6,1:6,1:nsta)/cmplx(0,omega) 
        stresssngl(1:6,1:6,1:nsta) = stress(1:6,1:6,1:nsta)
        displacementsngl(1:3,1:6,1:nsta) = displacement(1:3,1:6,1:nsta)

     endif

     coutfile = "out_stress_SH.bin"
     open(1,file=coutfile,status='unknown',form='unformatted',access='direct',recl=6*6*nsta*2*kind(0e0))
     write(1,rec=i+1) stresssngl(1:6,1:6,1:nsta)
     close(1)

     coutfile = "out_displ_SH.bin"
     open(1,file=coutfile,status='unknown',form='unformatted',access='direct',recl=3*6*nsta*2*kind(0e0))
     write(1,rec=i+1) displacementsngl(1:3,1:6,1:nsta)
     close(1)          

  enddo !!-- worker loop

  call date_and_time(datex,timex)
  write(*,'("end ",a4,"-",a2,"-",a2,"T",a2,":",a2,":",a4)') datex(1:4),datex(5:6),datex(7:8), &
        timex(1:2),timex(3:4),timex(5:8)

stop

end program myTraSH
