program myTraFFT
use hdf5
use mspline
use menv
use mrot
implicit none
character(40) :: nodesfile,normsfile,indexfile,freqfile,attr,f_text
integer :: i,j,k,imin,imax,i_source,i_channel,i_geocentric,i_filter,i_basis,imt,nsta,lsmooth,np0,np1,n_dcm,n_spl,n_t,n_tdwn,ntime,nidx,nsta_text
real(kind(0d0)) :: tlen,fmin,fmax,r0,lat0,lon0,samplingHz,d_t,norm
integer, dimension(:), allocatable :: updown,indx
real(kind(0d0)), dimension(3) :: normal_real
real(kind(0d0)), dimension(6) :: mt
real(kind(0d0)), dimension(:), allocatable :: stla,stlo,r,theta,phi,weight,bspl,ddiag
real(kind(0e0)), dimension(:,:), allocatable :: ygt
real(kind(0d0)), dimension(:,:), allocatable :: diag
real(kind(0d0)), dimension(:,:,:), allocatable :: fields
complex(kind(0d0)), dimension(:,:), allocatable :: normal,gt
complex(kind(0e0)), dimension(:,:,:), allocatable :: tmpsngl
complex(kind(0e0)), dimension(:,:,:,:), allocatable :: stresssngl,displacementsngl

call myinput(nodesfile,normsfile,indexfile,mirrorfile,tlen,imin,imax,fmin,fmax,r0,lat0,lon0,i_source,i_channel,i_geocentric,i_filter,i_basis,mt)

open(1,file=nodesfile)
read(1,*) nsta
allocate(r(1:nsta))
allocate(theta(1:nsta))
allocate(phi(1:nsta))
allocate(stla(1:nsta))
allocate(stlo(1:nsta))
allocate(updown(1:nsta))
do i = 1,nsta
  read(1,*) r(i),stla(i),stlo(i),updown(i)
  r(i) = 6371.d0-r(i)
  if (i_geocentric.eq.1) call translat(stla(i),stla(i))
  call calthetaphi(lat0,lon0,stla(i),stlo(i),theta(i),phi(i))
enddo
close(1)

open(1,file=normsfile)
allocate(weight(1:nsta),normal(1:nsta,1:3))
read(1,*) nsta
do i = 1,nsta
  read(1,*) weight(i),normal_real(1:3)
  norm = dsqrt(normal_real(1)**2+normal_real(2)**2+normal_real(3)**2)
  normal(i,1) = dcmplx(normal_real(1)/norm)
  normal(i,2) = dcmplx(normal_real(2)/norm)
  normal(i,3) = dcmplx(normal_real(3)/norm)
enddo
close(1)

open(1,file=indexfile)
read(1,*) nidx
allocate(indx(1:nidx))
do i = 1,nidx
  read(1,*) indx(i)
enddo
close(1)

!!! GB GB
!samplingHz = 20.d0
!np0 = imax
!call lsmoothfinder(tlen,np0,samplingHz,lsmooth)
!np1 = imax*lsmooth
!n_t = 2*np1
!samplingHz = dble(n_t)/tlen
samplingHz = 10.d0
call lsmoothfinder(tlen,imax,samplingHz,lsmooth)
i = 1
do while (i<lsmooth)
  i = i*2
enddo
lsmooth = i
np1 = 1
do while (np1<imax)
   np1 = np1*2
enddo
np1 = np1*lsmooth
n_t = 2*np1
samplingHz = dble(n_t)/tlen
write(*,*) "np1=",np1
write(*,*) "samplingHz=",samplingHz
write(*,*) "lsmooth=",lsmooth
!!! GB GB

allocate(tmpsngl(1:6,1:6,1:nsta))
allocate(stresssngl(1:6,1:6,1:nsta,imin:imax))
allocate(displacementsngl(1:3,1:6,1:nsta,imin:imax))

stresssngl=cmplx(0.e0)
displacementsngl=cmplx(0.e0)

do i = imin,imax
  ! PSV calculation
  freqfile = "out_stress_PSV.bin"
  open(1,file=freqfile,status='unknown',form='unformatted',access='direct',recl=6*6*nsta*2*kind(0e0))
  read(1,rec=i+1) tmpsngl(1:6,1:6,1:nsta)
  stresssngl(1:6,1:6,1:nsta,i) = stresssngl(1:6,1:6,1:nsta,i)+tmpsngl(1:6,1:6,1:nsta)
  close(1)   
  freqfile = "out_displ_PSV.bin"
  open(1,file=freqfile,status='unknown',form='unformatted',access='direct',recl=3*6*nsta*2*kind(0e0))
  read(1,rec=i+1) tmpsngl(1:3,1:6,1:nsta)
  displacementsngl(1:3,1:6,1:nsta,i) = displacementsngl(1:3,1:6,1:nsta,i)+tmpsngl(1:3,1:6,1:nsta)
  close(1)     
  ! SH calculation
  freqfile = "out_stress_SH.bin"
  open(1,file=freqfile,status='unknown',form='unformatted',access='direct',recl=6*6*nsta*2*kind(0e0))
  read(1,rec=i+1) tmpsngl(1:6,1:6,1:nsta)
  stresssngl(1:6,1:6,1:nsta,i) = stresssngl(1:6,1:6,1:nsta,i)+tmpsngl(1:6,1:6,1:nsta)
  close(1)
  freqfile = "out_displ_SH.bin"
  open(1,file=freqfile,status='unknown',form='unformatted',access='direct',recl=3*6*nsta*2*kind(0e0))
  read(1,rec=i+1) tmpsngl(1:3,1:6,1:nsta)
  displacementsngl(1:3,1:6,1:nsta,i) = displacementsngl(1:3,1:6,1:nsta,i)+tmpsngl(1:3,1:6,1:nsta)
  close(1)
enddo
deallocate(tmpsngl)

allocate(fields(1:6,0:n_t-1,1:nsta))
allocate(gt(1:9,0:n_t-1))
allocate(ygt(1:9,0:n_t-1))
do i = 1,nsta
  gt = cmplx(0.d0)
  do imt = 1,6
    if (i_source==1.and.imt>3) cycle
    gt(1:6,imin:imax) = gt(1:6,imin:imax)+stresssngl(1:6,imt,i,imin:imax)*mt(imt)
    gt(7:9,imin:imax) = gt(7:9,imin:imax)+displacementsngl(1:3,imt,i,imin:imax)*mt(imt)
  enddo
  call dsm_to_cartesian(r(i),stla(i),stlo(i),phi(i),imin,imax,n_t,normal(i,1:3),gt(1:9,0:n_t-1),i_basis)
  ygt = 0.e0
  call tensorFFT_real(9,imin,np1,gt,ygt,tlen,i_channel)
  fields(4:6,0:n_t-1,i) = dble(ygt(1:3,0:n_t-1))
  fields(1:3,0:n_t-1,i) = dble(ygt(7:9,0:n_t-1))
  if (i_filter==1) then
    do j = 1,3
      call bwfilt(dble(ygt(j,0:n_t-1)),fields(j+3,0:n_t-1,i),1.d0/samplingHz,n_t,1,2,fmin,fmax)
    enddo
    do j = 7,9
      call bwfilt(dble(ygt(j,0:n_t-1)),fields(j-6,0:n_t-1,i),1.d0/samplingHz,n_t,1,2,fmin,fmax)
    enddo
  endif
enddo
deallocate(gt,ygt)
deallocate(stresssngl,displacementsngl)

!!! NOBU !!!
!!! dump 10 first receivers in text format
nsta_text = min(10,nsta)
do i = 1,nsta_text

  write(f_text,'("out",i1.1,".txt")') i-1
  open(1,file=f_text)
  do j = 0,n_t-1
    write(1,'(f12.4,6e14.6)') j*1.d0/samplingHz,fields(7:9,j,i)
  enddo
  close(1)

enddo
!!! NOBU !!!

call create_mirror_file()

n_spl = 5
n_dcm = 10
d_t = n_dcm*1.d0/samplingHz
n_tdwn = int(n_t/n_dcm)+(n_spl+1)

attr = "Delta_t"
call create_attr_r(attr,d_t)
attr = "Num_t"
call create_attr_i(attr,n_tdwn)

allocate(bspl((n_spl+1)*n_dcm+1))
call bmn(bspl,n_dcm,n_spl)

allocate(buffer(6,nidx))
allocate(buffer_krn(6,nidx,n_spl+1))
buffer_krn = 0.d0

do ntime = 0,n_t-1
  if (n_dcm==1) then
    do k = 1,nidx
      buffer(:,k) = fields(:,ntime,indx(k))!*weight(indx(k))
    enddo
    call append_mirror_sl()
  else
    do i = 1,n_spl+1
      j = mod(ntime,n_dcm)+(n_spl+1-i)*n_dcm+1
      do k = 1,nidx
        buffer_krn(:,k,i) = buffer_krn(:,k,i)+bspl(j)*fields(:,ntime,indx(k))!*weight(indx(k))
      enddo
    enddo
    if (mod(ntime,n_dcm)==n_dcm-1) then
      buffer(:,:) = buffer_krn(:,:,1)
      call append_mirror_sl()
      do i = 1,n_spl
        buffer_krn(:,:,i) = buffer_krn(:,:,i+1)
      enddo
      buffer_krn(:,:,n_spl+1) = 0.d0
    endif
    if (ntime==n_t-1) then
      do i = 1,n_spl+1
        buffer(:,:) = buffer_krn(:,:,i)
        call append_mirror_sl()
      enddo
      allocate(diag(n_spl+1,n_tdwn),ddiag(n_tdwn))
      call bmdiag(n_dcm,n_spl,n_tdwn,n_t,diag)
      call bchfac(diag,n_spl+1,n_tdwn,ddiag)
      call bchslv_sl(diag,n_spl+1,n_tdwn)
      deallocate(diag,ddiag)
    endif
  endif
enddo

contains

subroutine create_mirror_file()
  implicit none
  character(len=40) :: dname
  integer(HID_T) :: h5_gzip_prop_2d,fid,dset_id
  integer :: hdferr

  call h5open_f(hdferr)
  call h5pcreate_f(H5P_DATASET_CREATE_F,h5_gzip_prop_2d,hdferr)
  call h5pset_deflate_f(h5_gzip_prop_2d,5,hdferr)
  call h5fcreate_f(mirrorfile,H5F_ACC_TRUNC_F,fid,hdferr)
  dname = "Fields_sl"
  call create_dset_2d(fid,trim(adjustl(dname)),H5T_IEEE_F64LE,int(6,HSIZE_T),int(H5S_UNLIMITED_F,HSIZE_T),dset_id)
  call h5dclose_f(dset_id,hdferr)
  dname = "Fields_fl"
  call create_dset_2d(fid,trim(adjustl(dname)),H5T_IEEE_F64LE,int(2,HSIZE_T),int(H5S_UNLIMITED_F,HSIZE_T),dset_id)
  call h5dclose_f(dset_id,hdferr)
  call h5fclose_f(fid,hdferr)

end subroutine create_mirror_file

subroutine create_attr_r(aname,attr)
  implicit none
  character(len=40), intent(in) :: aname
  real(kind(0d0)), intent(in) :: attr
  character(len=40) :: dname
  integer(HSIZE_T), dimension(1) :: dims,data_dims
  integer(HID_T) :: fid,dset_id,memspace,type_id,attr_id
  integer :: hdferr

  dims = (/1/)
  call h5fopen_f(mirrorfile,H5F_ACC_RDWR_F,fid,hdferr)
  call h5screate_f(H5S_SCALAR_F,memspace,hdferr)
  call h5acreate_f(fid,aname,H5T_NATIVE_DOUBLE,memspace,attr_id,hdferr)
  call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,attr,dims,hdferr)
  call h5aclose_f(attr_id,hdferr)
  call h5sclose_f(memspace,hdferr)
  call h5fclose_f(fid,hdferr)

end subroutine create_attr_r

subroutine create_attr_i(aname,attr)
  implicit none
  character(len=40), intent(in) :: aname
  integer, intent(in) :: attr
  character(len=40) :: dname
  integer(HSIZE_T), dimension(1) :: dims,data_dims
  integer(HID_T) :: fid,dset_id,memspace,type_id,attr_id
  integer :: hdferr

  dims = (/1/)
  call h5fopen_f(mirrorfile,H5F_ACC_RDWR_F,fid,hdferr)
  call h5screate_f(H5S_SCALAR_F,memspace,hdferr)
  call h5acreate_f(fid,aname,H5T_NATIVE_INTEGER,memspace,attr_id,hdferr)
  call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attr,dims,hdferr)
  call h5aclose_f(attr_id,hdferr)
  call h5sclose_f(memspace,hdferr)
  call h5fclose_f(fid,hdferr)

end subroutine create_attr_i

subroutine append_mirror_sl()
  implicit none
  character(len=40) :: dname
  integer(HID_T) :: fid,dset_id
  integer :: hdferr

  call h5fopen_f(mirrorfile,H5F_ACC_RDWR_F,fid,hdferr)
  dname = "Fields_sl"
  call h5dopen_f(fid,trim(adjustl(dname)),dset_id,hdferr)
  call append_dataset_2d(dset_id,buffer,hdferr)
  call h5dclose_f(dset_id,hdferr)
  call h5fclose_f(fid,hdferr)

end subroutine append_mirror_sl

subroutine append_dataset_2d(dset_id,arr,hdferr)
  implicit none
  integer(HID_T), intent(in) :: dset_id
  real(kind(0d0)), dimension(:,:), intent(in) :: arr
  integer, intent(out) :: hdferr
  integer(HSIZE_T), dimension(2) :: dims,maxdims,offset,dsize
  integer(HID_T) :: memspace,filespace

  dims(1) = size(arr,1)
  dims(2) = size(arr,2)
  call h5screate_simple_f(2,dims,memspace,hdferr)
  call h5dget_space_f(dset_id,filespace,hdferr)
  call h5sget_simple_extent_dims_f(filespace,dims,maxdims,hdferr)
  call h5sclose_f(filespace,hdferr)
  dsize(1) = size(arr,1)
  dsize(2) = dims(2)+size(arr,2)
  call h5dextend_f(dset_id,dsize,hdferr)
  call h5dget_space_f(dset_id,filespace,hdferr)
  offset(1) = 0
  offset(2) = dims(2)
  dims(1) = size(arr,1)
  dims(2) = size(arr,2)
  call H5Sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,hdferr)
  call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,arr,dims,hdferr,memspace,filespace)
  call H5Sclose_f(filespace,hdferr)
  call H5Sclose_f(memspace,hdferr)

end subroutine append_dataset_2d

subroutine create_dset_2d(f_id,dname,dtype,d1,d2,dset_id)
  character(len=*), intent(in) :: dname
  integer(HID_T), intent(in) :: f_id,dtype
  integer(HID_T), intent(out) :: dset_id
  integer(HSIZE_T), intent(in) :: d1,d2
  integer(HSIZE_T), dimension(2) :: dims,chunk,maxdims
  integer(HID_T) :: space_id,prop_id
  integer :: hdferr

  dims(1) = d1
  dims(2) = d2
  maxdims(1) = d1
  maxdims(2) = d2
  chunk(1) = min(d1,256*1024_HSIZE_T)
  if (d2==H5S_UNLIMITED_F) then
    chunk(2) = 64
    dims(2) = 0
  else
    chunk(2) = max(1_HSIZE_T,min(d2,int(256*1024/chunk(1),HSIZE_T)))
  endif
  call h5screate_simple_f(2,dims,space_id,hdferr,maxdims)
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdferr)
  if ((d1*d2).gt.128.or.d2==H5S_UNLIMITED_F) then
    call h5pset_chunk_f(prop_id,2,chunk,hdferr)
    if (dtype/=H5T_IEEE_F32LE.and.dtype/=H5T_IEEE_F64LE) then
      call h5pset_deflate_f(prop_id,5,hdferr)
      call h5pset_shuffle_f(prop_id,hdferr)
    endif
  end if
  call h5dcreate_f(f_id,dname,dtype,space_id,dset_id,hdferr,prop_id)
  call h5pclose_f(prop_id,hdferr)
  call h5sclose_f(space_id,hdferr)

end subroutine create_dset_2d

subroutine read_datasubset_2d(dset_id,offset,count,arr,hdferr)
  implicit none
  integer(HID_T), intent(in) :: dset_id
  integer(HSIZE_T), dimension(2), intent(in) ::  offset,count
  real(kind(0d0)), dimension(:,:), allocatable, intent(out) :: arr
  integer, intent(out) :: hdferr
  integer(HSIZE_T), dimension(2) :: stride,block,dimsm
  integer(HID_T) :: dataspace,memspace

  stride = (/1,1/)
  block = (/1,1/)
  dimsm = count
  call h5dget_space_f(dset_id,dataspace,hdferr)
  call h5sselect_hyperslab_f(dataspace,H5S_SELECT_SET_F,offset,count,hdferr,stride,block)
  call h5screate_simple_f(2,dimsm,memspace,hdferr)
  allocate(arr(dimsm(1),dimsm(2)))
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,arr,dimsm,hdferr,memspace,dataspace)
  call h5sclose_f(dataspace,hdferr)
  call h5sclose_f(memspace,hdferr)

end subroutine read_datasubset_2d

subroutine write_datasubset_2d(dset_id,offset,count,arr,hdferr)
  implicit none
  integer(HID_T), intent(in) :: dset_id
  integer(HSIZE_T), dimension(2), intent(in) :: offset,count
  real(kind(0d0)), dimension(:,:), intent(in) :: arr
  integer, intent(out) :: hdferr
  integer(HSIZE_T), dimension(2) :: stride,block,dimsm
  integer(HID_T) :: dataspace,memspace

  stride = (/1,1/)
  block = (/1,1/)
  dimsm = count
  call h5dget_space_f(dset_id,dataspace,hdferr)
  call h5sselect_hyperslab_f(dataspace,H5S_SELECT_SET_F,offset,count,hdferr,stride,block)
  call h5screate_simple_f(2,dimsm,memspace,hdferr)
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,arr,dimsm,hdferr,memspace,dataspace)
  call h5sclose_f(dataspace,hdferr)
  call h5sclose_f(memspace,hdferr)

end subroutine write_datasubset_2d

subroutine read_mirror_sl(n)
  implicit none
  integer, intent(in) :: n
  character (len=40) :: dname
  integer(HID_T) :: fid,dset_id
  integer :: hdferr
  integer(HSIZE_T), dimension(2) :: offset,dcount

  offset = (/0,n*size(buffer,2)/)
  dcount = (/6,size(buffer,2)/)
  call h5fopen_f(mirrorfile,H5F_ACC_RDONLY_F,fid,hdferr)
  dname = "Fields_sl"
  call h5dopen_f(fid,trim(dname),dset_id,hdferr)
  call read_datasubset_2d(dset_id,offset,dcount,buffer,hdferr)
  call h5dclose_f(dset_id,hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine read_mirror_sl

subroutine write_mirror_sl(n)
  implicit none
  integer, intent(in) :: n
  character (len=40) :: dname
  integer(HID_T) :: fid,dset_id
  integer :: hdferr
  integer(HSIZE_T), dimension(2) :: offset,dcount

  offset = (/0,n*size(buffer,2)/)
  dcount = (/6,size(buffer,2)/)
  call h5fopen_f(mirrorfile,H5F_ACC_RDWR_F,fid,hdferr)
  dname = "Fields_sl"
  call h5dopen_f(fid,trim(dname),dset_id,hdferr)
  call write_datasubset_2d(dset_id,offset,dcount,buffer,hdferr)
  call h5dclose_f(dset_id,hdferr)
  call h5fclose_f(fid,hdferr)

end subroutine write_mirror_sl

subroutine bchslv_sl(w,n,m)
  implicit none
  integer, intent(in) :: m,n
  real(kind(0d0)), intent(in) :: w(n,m)
  integer :: i,j

  !! Forward substitution, Solve L*Y = B.
  do i = 1,m-n+1
    if (i==1) then
      do j = 0,n-1
        call read_mirror_sl(j+i-1)
        buffer_krn(:,:,j+1) = buffer(:,:)
      enddo
    else
      do j = 0,n-2
        buffer_krn(:,:,j+1) = buffer_krn(:,:,j+2)
      enddo
      call read_mirror_sl(i+n-2)
      buffer_krn(:,:,n) = buffer(:,:)
    endif
    do j = 1,n-1
      buffer_krn(:,:,j+1) = buffer_krn(:,:,j+1)-w(j+1,i)*buffer_krn(:,:,1)
    enddo
    buffer(:,:) = buffer_krn(:,:,1)
    call write_mirror_sl(i-1)
  enddo
  do i = m-n+2,m
    do j = 0,n-2
      buffer_krn(:,:,j+1) = buffer_krn(:,:,j+2)
    enddo
    do j = 1,m-i
      buffer_krn(:,:,j+1) = buffer_krn(:,:,j+1)-w(j+1,i)*buffer_krn(:,:,1)
    enddo
    buffer(:,:) = buffer_krn(:,:,1)
    call write_mirror_sl(i-1)
  enddo
  !! Back substitution, Solve L'*X = D**(-1)*Y.
  do i = m,m-n+2,-1
    do j = m-i,1,-1
      buffer_krn(:,:,j+1) = buffer_krn(:,:,j)
    enddo
    call read_mirror_sl(i-1)
    buffer_krn(:,:,1) = buffer(:,:)*w(1,i)
    do j = 1,m-i
      buffer_krn(:,:,1) = buffer_krn(:,:,1)-w(j+1,i)*buffer_krn(:,:,j+1)
    enddo
    buffer(:,:) = buffer_krn(:,:,1)
    call write_mirror_sl(i-1)
  enddo
  do i = m-n+1,1,-1
    do j = n-1,1,-1
      buffer_krn(:,:,j+1) = buffer_krn(:,:,j)
    enddo
    call read_mirror_sl(i-1)
    buffer_krn(:,:,1) = buffer(:,:)*w(1,i)
    do j = 1,n-1
      buffer_krn(:,:,1) = buffer_krn(:,:,1)-w(j+1,i)*buffer_krn(:,:,j+1)
    enddo
    buffer(:,:) = buffer_krn(:,:,1)
    call write_mirror_sl(i-1)
  enddo
  return

end subroutine bchslv_sl

end program myTraFFT
