subroutine calbvecphi0_stock( l,thetadeg,plmDir,bvec,bvecdt,bvecdp)
  
  implicit none
  character(*) :: plmDir
  character(120) :: coutfile
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  integer  :: l,m,i,j
  real(kind(0d0)) :: theta,thetadeg,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

 
  write(coutfile, '(I8, ".","PLM")') int(thetadeg*100000.d0)
  do j = 1,8
     if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
  enddo
     
  coutfile = plmDir//"/"//coutfile

  open(1,file=coutfile,status='old',form='unformatted',access='direct', &
       recl=kind(0d0)*12)
  read(1,rec=l+1)plm(1,0),plm(1,1),plm(1,2),plm(1,3), &
       plm(2,0),plm(2,1),plm(2,2),plm(2,3), &
       plm(3,0),plm(3,1),plm(3,2),plm(3,3)
  close(1)

  theta =  thetadeg/180.d0*pi

  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     plmdt = dble(m) * x / sin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta) * coef * plm(1,m) 
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) = - coef * plmdt 
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     ! calculate derivatives
     bvecdt(1,m)  = dcmplx( 0.d0 )
     bvecdt(1,-m) = dcmplx( 0.d0 )
     bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt / dsin(theta) &
          - x / ( 1 - x * x ) * plm(1,m) ) * coef 
     bvecdt(2,-m) = dcmplx( conjg( bvecdt(2,m) ))
     bvecdt(3,m) = ( x / dsin(theta) * plmdt - dble(m) * dble(m)/(1-x*x) *plm(1,m) &
          &           + xl2 * plm(1,m) ) * coef 
     bvecdt(3,-m) = dcmplx(conjg( bvecdt(3,m)) )
     bvecdp(1,m)  = dcmplx( 0.d0 )
     bvecdp(1,-m) = dcmplx( 0.d0 )
     bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m) * coef 
     bvecdp(2,-m) = dcmplx(conjg( bvecdp(2,m)) )
     bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef 
     bvecdp(3,-m) = dcmplx(conjg( bvecdp(3,m)) )

     if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo

  return
end subroutine calbvecphi0_stock



subroutine calbvecphi0( l,theta,plm,bvec,bvecdt,bvecdp)
  
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  integer  :: l,m,i,j
  real(kind(0d0)) :: theta,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)
  do m=0,min0(l,3)
     call calplm( l,m,x,plm(1:3,m))
  enddo 

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     plmdt = dble(m) * x / sin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta) * coef * plm(1,m) 
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) = - coef * plmdt 
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     ! calculate derivatives
     bvecdt(1,m)  = dcmplx( 0.d0 )
     bvecdt(1,-m) = dcmplx( 0.d0 )
     bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt / dsin(theta) &
          - x / ( 1 - x * x ) * plm(1,m) ) * coef 
     bvecdt(2,-m) = dcmplx( conjg( bvecdt(2,m) ))
     bvecdt(3,m) = ( x / dsin(theta) * plmdt - dble(m) * dble(m)/(1-x*x) *plm(1,m) &
          &           + xl2 * plm(1,m) ) * coef 
     bvecdt(3,-m) = dcmplx(conjg( bvecdt(3,m)) )
     bvecdp(1,m)  = dcmplx( 0.d0 )
     bvecdp(1,-m) = dcmplx( 0.d0 )
     bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m) * coef 
     bvecdp(2,-m) = dcmplx(conjg( bvecdp(2,m)) )
     bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef 
     bvecdp(3,-m) = dcmplx(conjg( bvecdp(3,m)) )

     if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo



  return
end subroutine calbvecphi0



subroutine calplm( l,m,x,plm )
  implicit none
  integer :: l,m,i
  real(kind(0d0)) :: x,plm(1:3),pmm,somx2,fact

  if ((m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.d0)) stop 'bad arguments'
  if ( l.eq.m ) then
     pmm = 1.d0
     if ( m.gt.0 ) then
        somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
        fact = 1.d0
        do i=1,m
           pmm = -pmm * fact * somx2
           fact = fact + 2.d0
        enddo
     endif
     plm(3) = 0.d0
     plm(2) = 0.d0
     plm(1) = pmm
  else
     plm(3) = plm(2)
     plm(2) = plm(1)
     if ( l.eq.m+1 ) then
        plm(1) = x * dble(2*m+1) * plm(2)
     else
        !print *, l,m,x
        plm(1) = (x*dble(2*l-1) * plm(2)-dble(l+m-1) * plm(3) )/dble(l-m)
        !print *, plm(1)
     endif
  endif


end subroutine calplm


subroutine calbveczero( l,bvec )
  
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  
  integer  :: l,m,i
  real(kind(0d0)) :: fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  real(kind(0d0)) :: xl2


  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,1)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m)) *  xl2 *coef / 2.d0
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) =  -dcmplx(dble(m),0.d0) * xl2 * coef / 2.d0
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
     endif
  enddo



  return
end subroutine calbveczero



subroutine calbvec( l,theta,phi,plm,bvec,bvecdt,bvecdp )
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c Evaluating the value of toroidal harmonics (fully normalized)
  !c at each station whose latitude and longitude are theta and phi.
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  real(kind(0d0)),parameter :: pi=3.1415926535897932d0 
  
  integer ::  l,m,i
  real(kind(0d0)):: theta,phi,x,plm(3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(3,-2:2),expimp
  complex(kind(0d0)) :: bvecdt(3,-2:2),bvecdp(3,-2:2)
  real(kind(0d0)) :: plmdt,xl2
  
  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)
  do m=0,min0(l,3)
     call calplm( l,m,x,plm(1,m) )
  enddo
  do m=0,min0(l,2)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     expimp = cdexp( dcmplx( 0.d0, dble(m)*phi ) )
     plmdt = dble(m) * x / dsin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta ) * coef * plm(1,m) * expimp
     bvec(2,-m) = dconjg( bvec(2,m) )
     bvec(3,m) = - coef * plmdt * expimp
     bvec(3,-m) = dconjg( bvec(3,m) )
     ! calculate derivatives
     bvecdt(1,m)  = dcmplx( 0.d0 )
     bvecdt(1,-m) = dcmplx( 0.d0 )
     bvecdt(2,m)  = dcmplx( 0.d0, dble(m) )* ( plmdt / dsin(theta)- x/(1-x*x)*plm(1,m))*coef*expimp
     bvecdt(2,-m) = dconjg( bvecdt(2,m) )
     bvecdt(3,m) = (x/dsin(theta)*plmdt-dble(m)*dble(m)/(1-x*x)*plm(1,m)+xl2*plm(1,m))*coef*expimp
     bvecdt(3,-m) = dconjg( bvecdt(3,m) )
     bvecdp(1,m)  = dcmplx( 0.d0 )
     bvecdp(1,-m) = dcmplx( 0.d0 )
     bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m)*coef * expimp
     bvecdp(2,-m) = dconjg( bvecdp(2,m) )
     bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef * expimp
     bvecdp(3,-m) = dconjg( bvecdp(3,m) )
     
     if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo
  return
end subroutine calbvec

