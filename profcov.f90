module profcov
  use hdf5
  use precision
  implicit none
  integer, parameter, private :: wp=dp
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

  ! Hypsometric thickness of a layer defined by two pressure
  ! boundaries p1 and p2
  ! p1 - [Pa]
  ! p2 - [Pa]
  ! T  - [K]
  ! dz - [m]
  elemental function hypsometric(p1,p2,T) result(dz)
    real(wp), intent(in) :: p1,p2,T
    real(wp) :: dz
    real(wp), parameter :: &
         Rd=287.0_wp, & ! [J/K/kg]
         g0=9.81_wp     ! [m/s^2]

    dz=Rd*T/g0*(log(p1)-log(p2))
  end function hypsometric

  subroutine test_refraction(nu,p,temp,rho_H2O,theta,verbose,x)
    use bodhaine
    real, intent(in) :: p(:),temp(:)
    real(dp), intent(in) :: nu,rho_H2O
    real(dp), intent(inout) :: theta
    real(dp) :: dn,dn_prev
    integer :: i,m
    real(dp), parameter :: atm=101325d0,mbar=1d2
    real(dp) :: X_CO2,s,dz,z_tot
    real(dp), intent(out) :: x
    logical, intent(in) :: verbose

    m=size(p)

    ! 1 atm=101325 Pa
    ! 1 bar=1e5 Pa
    !nu=1e7_wp/11300
    !nu=1e7_wp/5000
    !rho_H2O=1e3_wp
    X_CO2=400.0_wp
    s=sin(theta)
    x=0
    dn_prev=0
    z_tot=0

    do i=1,m-1
       dn=n_air_minus_1(nu,p(i)*mbar/atm,real(temp(i),dp),rho_H2O,X_CO2)
       dz=hypsometric(p(i)*mbar,p(i+1)*mbar,real(temp(i)+temp(i+1),dp)/2)
       z_tot=z_tot+dz
       x=x+tan(s)*dz
       if (verbose) then
          write (*,'(I2,1X,F6.1,1X,F7.3,1X,EN22.12,3(1X,EN22.12))') i,p(i),temp(i),theta/degree,dn,dz,x
       end if
       s=(1+dn)/(1+dn_prev)*s
       theta=asin(s)
    end do
    if (verbose) then
       print *,z_tot
    end if
  end subroutine test_refraction

  subroutine test
    integer(hid_t) :: fid
    integer :: error
    integer(hsize_t) :: nlat,nlon,nlev,ntime,ntemp(4)
    real, allocatable :: lats(:),lons(:),levs(:),times(:),temp(:,:,:,:),mu_temp(:)
    real(dp), allocatable :: cov_temp(:,:)
    real(dp) :: x,theta,sx,sxx,ex,exx
    character(len=:), allocatable :: fn
    integer(hsize_t) :: ilat,ilat2,ilon,ilon2,ilev,itime,ilev1,ilev2,ncov,nx,ilambda,ih2o
    real(dp), parameter :: lambdas(*)=[5.0,15.0], & !,1.0,5.0,10.0,15.0]
                           h2os(*)=[1e2,5e4] !2e2,5e2,1e3,2e3,5e3,1e4,2e4,5e4]
    logical :: first

    first=.true.

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

    ! Refraction test
    ilat=(1+nlat)/2
    itime=1
    do ih2o=1,size(h2os)
    do ilambda=1,size(lambdas)
       sx=0
       sxx=0
       nx=0
       do ilat=(1+nlat)/2-90,(1+nlat)/2+90
          do ilon=1,nlon
             x=0
             theta=85*degree
             call test_refraction(1e4/lambdas(ilambda),levs,temp(ilon,ilat,:,itime),&
                  h2os(ih2o),theta,first,x)
             first=.false.
             sx=sx+x
             sxx=sxx+x*x
             nx=nx+1
             ! write (*,'(2(1X,EN22.12))') theta,x
          end do
       end do
       ex=sx/nx
       exx=sxx/nx
       write (*,'(F8.3,1X,ES14.7,1X,A8,1X,EN22.12,A8,1X,EN22.12)') &
            lambdas(ilambda),h2os(ih2o),'X:mu=',ex,'sigma=',sqrt(exx-ex*ex)
    end do
 end do

    return

    ! Covariance
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
                   ilat2=modulo(ilat+1,nlat)+1
                   ilon2=ilon !modulo(ilon-1,nlon)+1
                   cov_temp(ilev1,ilev2)=cov_temp(ilev1,ilev2)+ &
                        (temp(ilon,ilat,ilev1,itime)-mu_temp(ilev1))* &
                        (temp(ilon2,ilat2,ilev2,itime)-mu_temp(ilev2))
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
