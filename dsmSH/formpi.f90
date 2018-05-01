subroutine simplesplit(istart,iend,n,imin,imax)
!!! this routine returns the imin(i) and imax(i) (i=1,..,n) separates istart-iend into n parts
!!! each part contains imin(i)-imax(i) (i=1,..,n), iend-istart+1=inum=n*deltai+remainder,
!!! remainder is included in i th irange
!!! iend=n*deltai+remainder+istart-1, istart=iend-remainder-n*deltai+1
!!! n*deltai=iend-(start+remainder)+1
!!! remain :: istart,istart+1,..,(istart+remainder-1)
!!! the others :: i=(istart+iamari),..,imax -- (deltai*n)
!!! i(i) = istart+remainder+(i-1)*deltai,..,istart+remainder+n*deltai-1
  implicit none 
  integer, intent(in) :: istart,iend
  integer, intent(in) :: n ! the number of processors
  integer, dimension(*), intent(out) :: imin,imax
  integer :: remainder
  integer :: inum
  integer :: deltai
  integer :: i

  inum = iend-istart+1
  remainder = mod(inum,n)
  deltai = (inum-remainder)/n
  imin(1) = istart
  imax(1) = istart+remainder-1+deltai
  do i = 2,n
    imin(i) = istart+remainder+(i-1)*deltai
    imax(i) = istart+remainder+i*deltai-1
  enddo
  return
end subroutine simplesplit
      
subroutine trianglesplit(istart,iend,n,imin,imax)
!!! this routine returns the imin(i) and imax(i) (i=1,...,n) separates istart-iend into n parts
!!! each part contains imin(i)-imax(i) (i=1,..,n), iend-istart+1 = inum 
!!! Assume that the cpu time t for i th omega is a*i then we divide a*0.5*i**2 into n parts.
!!! return imin imax which satisfy above assumption
  implicit none 
  integer, intent(in) :: istart,iend
  integer, intent(in) :: n ! the number of processors
  integer, dimension(*), intent(out) :: imin,imax
  integer :: remainder
  integer :: inum
  integer :: deltai
  integer :: i
  integer, dimension(0:n) :: x
  real(kind(0d0)) :: s ! 0.5*iend**2/n
  real(kind(0d0)) :: p

  inum = iend-istart+1
  s = iend*iend/n
  x(0) = istart
  do i = 1,n
     x(i) = s+x(i-1)**2
     x(i) = x(i)**0.5
  enddo
  do i = 1,n
     imin(i) = x(i-1)
     imax(i) = x(i)-1
  enddo
  imax(n) = iend
  return
end subroutine trianglesplit
