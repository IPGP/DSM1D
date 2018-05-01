program merging
implicit none
integer :: imin,imax,n_mrg,i,j,n_nodes_all,n_nodes_stack
integer, dimension(:), allocatable :: n_nodes
character(len=80), dimension(:), allocatable :: f_nodes,f_disp_psv,f_stress_psv
complex(kind(0e0)), dimension(:,:,:), allocatable :: disp,stress,tmp_disp,tmp_stress

open(1,file="config.mrg")
read(1,*) imin,imax
read(1,*) n_mrg
allocate(f_nodes(n_mrg),f_disp_psv(n_mrg),f_stress_psv(n_mrg))
do i = 1,n_mrg
  read(1,*) f_nodes(i),f_disp_psv(i),f_stress_psv(i)
enddo
close(1)

allocate(n_nodes(n_mrg))
n_nodes_all = 0
do i = 1,n_mrg
  open(1,file=f_nodes(i))
  read(1,*) n_nodes(i)
  n_nodes_all = n_nodes_all+n_nodes(i)
  close(1)
enddo

allocate(disp(3,6,n_nodes_all),stress(6,6,n_nodes_all))
open(2,file="out_displ_PSV.bin",form='unformatted',access='direct',recl=3*6*n_nodes_all*2*kind(0e0))
open(3,file="out_stress_PSV.bin",form='unformatted',access='direct',recl=6*6*n_nodes_all*2*kind(0e0))
do i = imin,imax
  write(*,*) i,imin,imax
  n_nodes_stack = 0
  do j = 1,n_mrg
    allocate(tmp_disp(3,6,n_nodes(j)),tmp_stress(6,6,n_nodes(j)))
    open(1,file=f_disp_psv(j),form='unformatted',access='direct',recl=3*6*n_nodes(j)*2*kind(0e0))
    read(1,rec=i+1) tmp_disp(1:3,1:6,1:n_nodes(j))
    disp(1:3,1:6,(n_nodes_stack+1):(n_nodes_stack+n_nodes(j))) = tmp_disp(1:3,1:6,1:n_nodes(j))
    close(1)
    open(1,file=f_stress_psv(j),form='unformatted',access='direct',recl=6*6*n_nodes(j)*2*kind(0e0))
    read(1,rec=i+1) tmp_stress(1:6,1:6,1:n_nodes(j))
    stress(1:6,1:6,(n_nodes_stack+1):(n_nodes_stack+n_nodes(j))) = tmp_stress(1:6,1:6,1:n_nodes(j))
    close(1)
    deallocate(tmp_disp,tmp_stress)
    n_nodes_stack = n_nodes_stack+n_nodes(j)
  enddo
  write(2,rec=i+1) disp(1:3,1:6,1:n_nodes_all)
  write(3,rec=i+1) stress(1:6,1:6,1:n_nodes_all)
enddo
close(2)
close(3)

end program merging
