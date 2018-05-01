module  mrot

public :: dsm_to_cartesian
private

contains

subroutine dsm_to_cartesian(r,lat,lon,phi,imin,imax,n,normal_xyz,gt,i_basis)
  implicit none
  integer, intent(in) :: imin,imax,n,i_basis
  real(kind(0d0)),intent(in) :: r,lat,lon,phi
  complex(kind(0d0)), dimension(1:3), intent(in) :: normal_xyz
  complex(kind(0d0)), dimension(1:9,0:n-1), intent(inout) :: gt
  integer :: i
  complex(kind(0d0)), dimension(1:3) :: normal_zne,normal_zrt
  complex(kind(0d0)), dimension(1:3) :: displ_zrt,displ_zne,displ_xyz
  complex(kind(0d0)), dimension(1:3) :: tract_zrt,tract_zne,tract_xyz
  complex(kind(0d0)), dimension(1:3,1:3) :: R_xyz_to_zne,R_zne_to_zrt
  complex(kind(0d0)), dimension(1:3,1:3) :: stress_zrt

  call matrix_xyz_to_zne(lat,lon,R_xyz_to_zne)
  call matrix_zne_to_zrt(phi,R_zne_to_zrt)
  call matrix_product(R_xyz_to_zne,normal_xyz,normal_zne)
  call matrix_product(R_zne_to_zrt,normal_zne,normal_zrt)

  do i = imin,imax
    stress_zrt(1,1) = gt(1,i)
    stress_zrt(1,2) = gt(4,i)
    stress_zrt(1,3) = gt(5,i)
    stress_zrt(2,1) = gt(4,i)
    stress_zrt(2,2) = gt(2,i)
    stress_zrt(2,3) = gt(6,i)
    stress_zrt(3,1) = gt(5,i)
    stress_zrt(3,2) = gt(6,i)
    stress_zrt(3,3) = gt(3,i)
    displ_zrt(1) = gt(7,i)
    displ_zrt(2) = gt(8,i)
    displ_zrt(3) = gt(9,i)
    call matrix_product(stress_zrt,normal_zrt,tract_zrt)
    if (i_basis==0) then
      gt(1,i) = tract_zrt(1)
      gt(2,i) = tract_zrt(2)
      gt(3,i) = tract_zrt(3)
    else
      call matrix_product_transp(R_zne_to_zrt,tract_zrt,tract_zne)
      call matrix_product_transp(R_zne_to_zrt,displ_zrt,displ_zne)
      if (i_basis==1) then
        gt(1,i) = tract_zne(1)
        gt(2,i) = tract_zne(2)
        gt(3,i) = tract_zne(3)
        gt(7,i) = displ_zne(1)
        gt(8,i) = displ_zne(2)
        gt(9,i) = displ_zne(3)
        if (i_basis==2) then
          call matrix_product_transp(R_xyz_to_zne,tract_zne,tract_xyz)
          call matrix_product_transp(R_xyz_to_zne,displ_zne,displ_xyz)
          gt(1,i) = tract_xyz(1)
          gt(2,i) = tract_xyz(2)
          gt(3,i) = tract_xyz(3)
          gt(7,i) = displ_xyz(1)
          gt(8,i) = displ_xyz(2)
          gt(9,i) = displ_xyz(3)
        endif
      endif
    endif
  enddo

end subroutine dsm_to_cartesian

subroutine matrix_xyz_to_zne(lat,lon,Rot)
  implicit none
  real(kind(0d0)), intent(in) :: lat,lon
  complex(kind(0d0)), dimension(1:3,1:3), intent(out) :: Rot
  real(kind(0d0)) :: dg2rd

  dg2rd = 3.1415926535897932d0/180.d0

  Rot(1,1) = dcmplx(dcos(lat*dg2rd)*dcos(lon*dg2rd))
  Rot(1,2) = dcmplx(dcos(lat*dg2rd)*dsin(lon*dg2rd))
  Rot(1,3) = dcmplx(dsin(lat*dg2rd))

  Rot(2,1) = dcmplx(-dsin(lat*dg2rd)*dcos(lon*dg2rd))
  Rot(2,2) = dcmplx(-dsin(lat*dg2rd)*dsin(lon*dg2rd))
  Rot(2,3) = dcmplx(dcos(lat*dg2rd))

  Rot(3,1) = dcmplx(-dsin(lon*dg2rd))
  Rot(3,2) = dcmplx(dcos(lon*dg2rd))
  Rot(3,3) = dcmplx(0.d0)

end subroutine matrix_xyz_to_zne

subroutine matrix_zne_to_zrt(phi,Rot)
  implicit none
  real(kind(0d0)), intent(in) :: phi
  complex(kind(0d0)), dimension(1:3,1:3), intent(out) :: Rot
  real(kind(0d0)) :: dg2rd

  dg2rd = 3.1415926535897932d0/180.d0

  Rot(1,1) = dcmplx(1.d0)
  Rot(1,2) = dcmplx(0.d0)
  Rot(1,3) = dcmplx(0.d0)

  Rot(2,1) = dcmplx(0.d0)
  Rot(2,2) = dcmplx(dcos(phi*dg2rd))
  Rot(2,3) = dcmplx(dsin(phi*dg2rd))

  Rot(3,1) = dcmplx(0.d0)
  Rot(3,2) = dcmplx(-dsin(phi*dg2rd))
  Rot(3,3) = dcmplx(dcos(phi*dg2rd))

end subroutine matrix_zne_to_zrt

subroutine matrix_product(Rot,x_in,x_out)
  implicit none
  complex(kind(0d0)), dimension(1:3), intent(in) :: x_in
  complex(kind(0d0)), dimension(1:3,1:3), intent(in) :: Rot(3,3)
  complex(kind(0d0)), dimension(1:3), intent(out) :: x_out
  integer :: i,j

  x_out(:) = 0.d0
  do i = 1,3
    do j = 1,3
      x_out(i) = x_out(i)+Rot(i,j)*x_in(j)
    enddo
  enddo

end subroutine matrix_product

subroutine matrix_product_transp(Rot,x_in,x_out)
  implicit none
  complex(kind(0d0)), dimension(1:3), intent(in) :: x_in
  complex(kind(0d0)), dimension(1:3,1:3), intent(in) :: Rot(3,3)
  complex(kind(0d0)), dimension(1:3), intent(out) :: x_out
  integer :: i,j

  x_out(:) = 0.d0
  do i = 1,3
    do j = 1,3
      x_out(i) = x_out(i)+Rot(j,i)*x_in(j)
    enddo
  enddo

end subroutine matrix_product_transp

end module mrot
