subroutine clPLM(plm,lmax,theta,theta_n)
  implicit none
 
  integer :: itheta, l, theta_n, lmax,m
  real(kind(0d0)) :: plm(1:3,0:3,1:theta_n,0:lmax),theta(1:theta_n)
  real(kind(0d0)) :: tmpthetainrad
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)) :: x, plmtmp(1:3,0:3)
  
  plm = 0.d0

  do itheta = 1, theta_n
     tmpthetainrad = theta(itheta)/180.d0*pi
     x = cos(tmpthetainrad)
     plmtmp = 0.d0
     do l = 0, lmax
        do m = 0, min0(l,3)
           call calplm(l,m,x,plmtmp(1:3,m))
        enddo
        plm(1:3,0:3,itheta,l) = plmtmp(1:3,0:3)
        !print *, plmtmp(1:3,0:3)
     enddo
  enddo
  return
end subroutine clPLM



subroutine clPLM_modified(l,theta,plmtmp,plm)
  implicit none
  integer :: l, m
  real(kind(0d0)) :: plmtmp(1:3,0:3), plm(1:3,0:3)
  real(kind(0d0)) :: x, theta

  x = cos(theta)
  do m = 0, min0(l,3)
    call calplm(l,m,x,plmtmp(1:3,m))
  enddo
  plm(1:3,0:3) = plmtmp(1:3,0:3)
  return
end subroutine clPLM_modified



subroutine caldvecphi0( l,theta,plm,bvec,bvecdt,bvecdp)
  
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  integer  :: l,m,i,j
  real(kind(0d0)) :: theta,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)
  !do m=0,min0(l,3)
  !   call calplm( l,m,x,plm(1:3,m))
  !enddo 

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     plmdt = dble(m) * x / sin( theta ) * plm(1,m) + plm(1,m+1)
     
     bvec(1,m)  = coef * plm(1,m) 
     bvec(1,-m) = dconjg( bvec(1,m) )
     bvec(2,m) = coef * plmdt 
     bvec(2,-m) = dconjg( bvec(2,m) )
     bvec(3,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta ) * coef * plm(1,m) 
     bvec(3,-m) = dconjg( bvec(3,m) )
     
     
     ! calculate derivatives
     
     bvecdt(1,m) = plmdt * coef 
     bvecdt(1,m) = dconjg( bvecdt(1,m) )
     bvecdt(2,m) = ( - x / dsin(theta) * plmdt + dble(m) * dble(m) / ( 1 - x * x ) * plm(1,m) - xl2 * plm(1,m) ) * coef
     bvecdt(2,-m) = dconjg( bvecdt(2,m) )
     bvecdt(3,m) = dcmplx( 0.d0, dble(m) ) * ( - x / ( 1- x * x ) * plm(1,m)+ 1.d0 / dsin(theta) * plmdt ) * coef 
     bvecdt(3,-m) = dconjg( bvecdt(3,m) )
     bvecdp(1,m) = dcmplx( 0.d0, dble(m) ) * plm(1,m) * coef 
     bvecdp(1,-m) = dconjg( bvecdp(1,m) )
     bvecdp(2,m) = dcmplx( 0.d0,dble(m) ) * plmdt * coef 
     bvecdp(2,-m) = dconjg( bvecdp(2,m) )
     bvecdp(3,m) = - dble(m) * dble(m) / dsin(theta)*plm(1,m)*coef
     bvecdp(3,-m) = dconjg( bvecdp(3,m) )
     if ( mod(m,2).eq.1 ) then
        bvec(1,-m) = - bvec(1,-m)
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(1,-m) = - bvecdt(1,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(1,-m) = - bvecdp(1,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
     
  enddo
  return
end subroutine caldvecphi0


subroutine caldvecphi0_good( l,theta,plm,bvec,bvecdt,bvecdp)
  
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
     
     bvec(1,m)  = coef * plm(1,m) 
     bvec(1,-m) = dconjg( bvec(1,m) )
     bvec(2,m) = coef * plmdt 
     bvec(2,-m) = dconjg( bvec(2,m) )
     bvec(3,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta ) * coef * plm(1,m) 
     bvec(3,-m) = dconjg( bvec(3,m) )
     
     
     ! calculate derivatives
     
     bvecdt(1,m) = plmdt * coef 
     bvecdt(1,m) = dconjg( bvecdt(1,m) )
     bvecdt(2,m) = ( - x / dsin(theta) * plmdt + dble(m) * dble(m) / ( 1 - x * x ) * plm(1,m) - xl2 * plm(1,m) ) * coef
     bvecdt(2,-m) = dconjg( bvecdt(2,m) )
     bvecdt(3,m) = dcmplx( 0.d0, dble(m) ) * ( - x / ( 1- x * x ) * plm(1,m)+ 1.d0 / dsin(theta) * plmdt ) * coef 
     bvecdt(3,-m) = dconjg( bvecdt(3,m) )
     bvecdp(1,m) = dcmplx( 0.d0, dble(m) ) * plm(1,m) * coef 
     bvecdp(1,-m) = dconjg( bvecdp(1,m) )
     bvecdp(2,m) = dcmplx( 0.d0,dble(m) ) * plmdt * coef 
     bvecdp(2,-m) = dconjg( bvecdp(2,m) )
     bvecdp(3,m) = - dble(m) * dble(m) / dsin(theta)*plm(1,m)*coef
     bvecdp(3,-m) = dconjg( bvecdp(3,m) )
     if ( mod(m,2).eq.1 ) then
        bvec(1,-m) = - bvec(1,-m)
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(1,-m) = - bvecdt(1,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(1,-m) = - bvecdp(1,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
     
  enddo
  return
end subroutine caldvecphi0_good




subroutine caldvec_dejaplm( l,theta,phi,plm,bvec,bvecdt,bvecdp )

  implicit none
  real(kind(0d0)),parameter :: pi=3.1415926535897932d0 
  integer :: l,m,i
  real(kind(0d0)):: theta,phi,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2),expimp
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)):: plmdt,xl2

  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)

  !do 100 m=0,min0(l,3)
  !call calplm( l,m,x,plm(1,m) )
  !enddo
  do m=0,min0(l,2)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     expimp = zexp( dcmplx( 0.d0, dble(m)*phi ) )
     plmdt = dble(m) * x / dsin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = coef * plm(1,m) * expimp
     bvec(1,-m) = dconjg( bvec(1,m) )
     bvec(2,m) = coef * plmdt * expimp
     bvec(2,-m) = dconjg( bvec(2,m) )
     bvec(3,m)  = dcmplx( 0.d0, dble(m))/dsin( theta )*coef*plm(1,m)*expimp
     bvec(3,-m) = dconjg( bvec(3,m))
     
     
     ! calculate derivatives
     bvecdt(1,m) = plmdt * coef * expimp
     bvecdt(1,-m) = dconjg( bvecdt(1,m) )
     bvecdt(2,m) =(-x/dsin(theta) * plmdt+dble(m)*dble(m)/(1-x*x)*plm(1,m)-xl2*plm(1,m))*coef*expimp
     bvecdt(2,-m) = dconjg( bvecdt(2,m) )
     bvecdt(3,m) = dcmplx( 0.d0, dble(m) ) * ( - x / ( 1- x * x ) * plm(1,m) +1.d0/dsin(theta)*plmdt)*coef*expimp
     bvecdt(3,-m) = dconjg( bvecdt(3,m) )
     bvecdp(1,m) = dcmplx( 0.d0, dble(m) ) * plm(1,m) * coef * expimp
     bvecdp(1,-m) = dconjg( bvecdp(1,m) )
     bvecdp(2,m) = dcmplx( 0.d0,dble(m) ) * plmdt * coef * expimp
     bvecdp(2,-m) = dconjg( bvecdp(2,m) )
     bvecdp(3,m) = - dble(m) * dble(m) / dsin(theta)*plm(1,m) * coef * expimp
     bvecdp(3,-m) = dconjg( bvecdp(3,m) )
     
     if ( mod(m,2).eq.1 ) then
        bvec(1,-m) = - bvec(1,-m)
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(1,-m) = - bvecdt(1,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(1,-m) = - bvecdp(1,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo
  return
end subroutine caldvec_dejaplm


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
        plm(1) = (x*dble(2*l-1) * plm(2)-dble(l+m-1) * plm(3) )/dble(l-m)
     endif
  endif

end subroutine calplm


subroutine caldveczero( l,bvec )
  
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  
  integer  :: l,m,i
  real(kind(0d0)) :: fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  real(kind(0d0)) :: xl2


  bvec = cmplx(0.d0)
  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,1)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )

     if(m.eq.0) then
        bvec(1,m) = cmplx(coef) 
     endif
     if(m.eq.1) then
        bvec(2,m) = dcmplx(dble(m),0.d0) * xl2 * coef / 2.d0
        bvec(2,-m) = conjg(bvec(2,m))
        bvec(3,m) = dcmplx( 0.d0, dble(m)) *  xl2 *coef / 2.d0
        bvec(3,-m) = conjg(bvec(3,m))
     endif

     if(mod(m,2).eq.1) then
        bvec(2,-m) = -bvec(2,-m)
        bvec(3,-m) = -bvec(3,-m)
     endif
  enddo
end subroutine caldveczero
