module profcov
  use hdf5
  implicit none
  integer, parameter :: dp=selected_real_kind(12), &
                        wp=dp
  real(wp), parameter :: &
       pi=3.1415926535897931159979634685441851615906_wp, &
       degree=pi/180.0_wp
  complex(wp), parameter :: imag=(0.0_wp,1.0_wp)
contains
  subroutine getarg(k,u)
    integer, intent(in) :: k
    integer :: m
    character(len=:), allocatable, intent(inout) :: u

    if (allocated(u)) deallocate(u)
    call get_command_argument(k,length=m)
    allocate(character(len=m) :: u)
    call get_command_argument(k,u)
  end subroutine getarg

  subroutine errchk(error)
    integer, intent(in) :: error
    if (error /= 0) error stop "Error check failed"
  end subroutine errchk

  subroutine h5readvec(fid,name,nx,xs,error)
    integer(hid_t), intent(in) :: fid
    integer(hid_t) :: did,sid,mid
    character(len=*), intent(in) :: name
    real, allocatable, intent(out) :: xs(:)
    integer(hsize_t), intent(out) :: nx
    integer :: rank
    integer(hsize_t) :: dims(1),maxdims(1)
    integer, intent(inout) :: error

    print *,'READ:',name
    call h5dopen_f(fid,name,did,error)
    call h5dget_space_f(did,sid,error)
    call h5sget_simple_extent_ndims_f(sid,rank,error)
    if (rank/=1) error stop 'Unexpected rank'
    call h5sget_simple_extent_dims_f(sid,dims,maxdims,error)
    call h5screate_simple_f(rank,dims,mid,error)
    nx=dims(1)
    allocate(xs(nx))
    call h5dread_f(did,H5T_NATIVE_REAL,xs,dims,error,mid,sid)
    call h5sclose_f(sid,error)
    call h5sclose_f(mid,error)
  end subroutine h5readvec

  subroutine h5reada4(fid,name,nx,xs,error)
    integer(hid_t), intent(in) :: fid
    integer(hid_t) :: did,sid,mid,aid
    character(len=*), intent(in) :: name
    real, allocatable, intent(out) :: xs(:,:,:,:)
    integer(hsize_t), intent(out) :: nx(4)
    integer :: rank
    integer(hsize_t) :: dims(4),maxdims(4),adim(1)
    real :: scal(1),off(1)
    integer, intent(inout) :: error

    print *,'READ:',name
    call h5dopen_f(fid,name,did,error)
    call h5dget_space_f(did,sid,error)
    call h5sget_simple_extent_ndims_f(sid,rank,error)
    if (rank/=4) error stop 'Unexpected rank'
    call h5sget_simple_extent_dims_f(sid,dims,maxdims,error)
    call h5screate_simple_f(rank,dims,mid,error)
    nx=dims
    allocate(xs(nx(1),nx(2),nx(3),nx(4)))
    call h5dread_f(did,H5T_NATIVE_REAL,xs,dims,error,mid,sid)

    adim(1)=1
    call h5aopen_by_name_f(did,'.','scale_factor',aid,error)
    call h5aread_f(aid,H5T_NATIVE_REAL,scal,adim,error)
    call h5aclose_f(aid,error)

    call h5aopen_by_name_f(did,'.','add_offset',aid,error)
    call h5aread_f(aid,H5T_NATIVE_REAL,off,adim,error)
    call h5aclose_f(aid,error)

    xs=xs*scal(1)+off(1)
    call h5sclose_f(sid,error)
    call h5sclose_f(mid,error)
  end subroutine h5reada4

  subroutine test
    integer(hid_t) :: fid
    integer :: error
    integer(hsize_t) :: nlat,nlon,nlev,ntime,ntemp(4)
    real, allocatable :: lats(:),lons(:),levs(:),times(:),temp(:,:,:,:),mu_temp(:)
    real(dp), allocatable :: cov_temp(:,:)
    character(len=:), allocatable :: fn
    integer(hsize_t) :: ilat,ilon,ilev,itime,ilev1,ilev2,ncov

    call getarg(1,fn)

    call h5open_f(error)
    call h5fopen_f(fn,H5F_ACC_RDONLY_F,fid,error)

    call h5readvec(fid,'latitude',nlat,lats,error)
    call h5readvec(fid,'longitude',nlon,lons,error)
    call h5readvec(fid,'level',nlev,levs,error)
    call h5readvec(fid,'time',ntime,times,error)
    write (*,'(4(1X,I8))') nlat,nlon,nlev,ntime
    call h5reada4(fid,'t',ntemp,temp,error)
    call h5fclose_f(fid,error)
    print *,'t:',shape(temp)
    print *,'t:',ntemp

    write (*,'(A2,1X,A8,1X,A8)') '#','p [mbar]','T [K]'
    allocate(mu_temp(nlev),cov_temp(nlev,nlev))
    do ilev=1,nlev
       mu_temp(ilev)=sum(temp(:,:,ilev,:))/real(nlat*nlon*ntime)
       write (*,'(I2,1X,F8.1,F8.3)') ilev,levs(ilev),mu_temp(ilev)
    end do
    cov_temp=0
    ! longitude,latitude,level,time
    do itime=1,ntime
       do ilat=1,nlat
          do ilon=1,nlon
             do ilev1=1,nlev
                do ilev2=ilev1,nlev
                   cov_temp(ilev1,ilev2)=cov_temp(ilev1,ilev2)+ &
                        (temp(ilon,ilat,ilev1,itime)-mu_temp(ilev1))* &
                        (temp(ilon,ilat,ilev2,itime)-mu_temp(ilev2))
                end do
             end do
          end do
       end do
    end do
    ncov=nlat*nlon*ntime
    cov_temp=cov_temp/ncov
    do ilev1=1,nlev
       write (*,'(I2,1X,ES15.7)') ilev1,sqrt(cov_temp(ilev1,ilev1))
       cov_temp(ilev1+1:nlev,ilev1)=cov_temp(ilev1,ilev1+1:nlev)
       do ilev2=1,nlev
          write (*,'(I2,1X,I2,1X,ES15.7)') &
               ilev1,ilev2, &
               cov_temp(ilev1,ilev2)/(sqrt(cov_temp(ilev1,ilev1))*sqrt(cov_temp(ilev2,ilev2)))
       end do
       write (*,'()')
    end do
  end subroutine test
end module profcov

program profcov_p
  use profcov
  implicit none

  call test
end program profcov_p
