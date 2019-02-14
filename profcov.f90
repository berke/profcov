module profcov
  use hdf5
  use precision
  implicit none
  integer, parameter, private :: wp=dp
  real(wp), parameter :: &
       pi=3.1415926535897931159979634685441851615906_wp, &
       degree=pi/180.0_wp
  complex(wp), parameter :: imag=(0.0_wp,1.0_wp)

  type, abstract :: boundary
   contains
     procedure(next_boundary_intf), deferred :: next
  end type boundary

  type, abstract :: interaction
     ! k0 : current element
     ! p0 : position within k0
     ! u0 : direction
     ! k1 : next element
     ! p1 : position at interface
     ! v1 : normal at interface point
     ! u1 : new direction
   contains
     procedure(interact_intf), deferred :: interact
  end type interaction

  abstract interface
       ! k0       : current element
       ! p0       : current position
       ! u0       : current direction normal
       ! k1       : next element
       ! p1       : position at next boundary following u0
       ! v1       : boundary normal at p1
       subroutine next_boundary_intf(this,k0,p0,u0,k1,p1,v1)
         import
         class(boundary) :: this
         integer, intent(in) :: k0
         real(dp), intent(in) :: p0(3),u0(3)
         integer, intent(out) :: k1
         real(dp), intent(out) :: p1(3),v1(3)
       end subroutine next_boundary_intf

     subroutine interact_intf(this,k0,u0,k1,v1,u1)
       import
       class(interaction) :: this
       integer, intent(in) :: k0,k1
       real(dp), intent(in) :: u0(3),v1(3)
       real(dp), intent(out) :: u1(3)
     end subroutine interact_intf
  end interface

  type, extends(interaction) :: snell_descartes
     real(dp), allocatable :: indices(:)
   contains
     procedure :: interact=>snde_interact
  end type snell_descartes

  type, extends(boundary) :: cylindrical
     real(dp), allocatable :: rs(:) ! Radiuses, decreasing
   contains
     procedure :: next=>cyl_next
  end type cylindrical
contains
  subroutine quadr(a,b,c,x,n)
    real(dp), intent(in) :: a,b,c
    real(dp), intent(out) :: x(2)
    integer, intent(out) :: n
    real(dp) :: delta,t

    !     ax^2+bx+c=0
    !     x^2 + xb/a + c/a=0             DIVIDE BY a
    ! (1) x^2 + xb/a = -c/a              SUBSTRACT c/a
    !
    ! (x+b/(2a))^2
    !   = x^2+b^2/(4a^2)+2xb/(2a)
    !   = x^2+b^2/(4a^2)+xb/a       SIMPLIFY
    !
    ! THUS
    !
    ! (2) (x+b/(2a))^2-b^2/(4a^2) = x^2+xb/a
    !
    ! SUBSTITUTING (2) INTO (1)
    !
    ! (3) (x+b/(2a))^2-b^2/(4a^2) = -c/a
    !     (x+b/(2a))^2 = b^2/(4a^2) - c/a
    !
    ! THEREFORE
    !
    !     (x+b/(2a)) = ±sqrt(b^2/(4a^2) - c/a)
    !              x = ±sqrt(-c/a+b^2/(4a^2))-b/(2a)
    !
    ! b^2/(4a^2)-c/a=
    ! b^2/(4a^2)-ac/aa=
    ! (b^2-ac)/(4a^2)
    delta=b*b/(4*a*a)-c/a
    if (delta<0) then
       ! No real solutions
       n=0
    else
       t=-b/(2*a)
       if (abs(delta)<epsilon(delta)) then
          ! One solution
          x(1)=t
          n=1
       else
          delta=sqrt(delta)
          x(1)=-delta+t
          x(2)=delta+t
          n=2
       end if
    end if
  end subroutine quadr

  subroutine cyl_next(this,k0,p0,u0,k1,p1,v1)
    class(cylindrical) :: this
    integer, intent(in) :: k0
    real(dp), intent(in) :: p0(3),u0(3)
    integer, intent(out) :: k1
    real(dp), intent(out) :: p1(3),v1(3)
    real(dp) :: r0,r1,ur,p0u0,lams(2)
    integer :: nlam

    ! Cylinder axis: z-axis
    ur=hypot(u0(1),u0(2))
    if (ur<=epsilon(ur)) error stop 'Motion parallel to axis'
    r0=hypot(p0(1),p0(2))
    p0u0=dot_product(p0(1:2),u0(1:2))
    if (p0u0<0) then
       ! Radius is decreasing
       if (k0 < size(this%rs)) then
          k1=k0+1
       else
          k1=k0
       end if
    else
       if (k0 > 1) then
          k1=k0-1
       else
          k1=k0
       end if
    end if
    ! TODO: Fix logic

    r1=this%rs(k1)
    ! ||p0(1:2)+lam*u0(1:2)||^2=r1^2
    ! (p0(1)+lam*u0(1))^2 + (p0(2)+lam*u0(2))^2 = r1^2
    !
    ! p0(1)^2 + lam^2*u0(1))^2 + 2*p0(1)*lam*u0(1)
    !   + p0(2)^2 + lam^2*u0(2))^2 + 2*p0(2)*lam*u0(2) = r1^2
    !
    ! (p0(1)^2 + p0(2)^2)
    !   + 2*(p0(1)*u0(1)+p0(2)*u0(2))*lam
    !   + lam^2*(u0(1)^2 + u0(2)^2) = r1^2
    !
    ! r0^2
    !   + 2*(p0(1)*u0(1)+p0(2)*u0(2))*lam
    !   + lam^2*(u0(1)^2 + u0(2)^2) = r1^2
    !
    !     2*(p0(1)*u0(1)+p0(2)*u0(2))*lam
    !   + lam^2*(u0(1)^2 + u0(2)^2) = r1^2-r0^2
    !
    ! ||p+Lu||=R
    ! <p+Lu|p+Lu>=R^2
    ! <p|p>+2<p|Lu>+<Lu|Lu>=R^2
    ! <p|p>+L 2<p|u>+L^2 <u|u>=R^2
    ! <p|p>+L 2<p|u>+L^2 <u|u>-R^2=0
    ! (<p|p>-R^2)+L 2<p|u>+L^2 <u|u>=0
    call quadr(ur*ur,2*p0u0,r0*r0-r1*r1,lams,nlam)
    if (nlam==0) then
       k1=-1
       return
    else
       p1=p0+lams(1)*u0
       if (dot_product(p1-p0,u0)<0) then
          if (nlam==1) then
             k1=-1
             return
          else
             p1=p0+lams(2)*u0
          end if
       end if
    end if

    v1(1)=p1(2)
    v1(2)=-p1(1) ! Orientation?
    v1(3)=0
    v1=v1/norm2(v1)
  end subroutine cyl_next

  subroutine snde_interact(this,k0,u0,k1,v1,u1)
    class(snell_descartes) :: this
    integer, intent(in) :: k0,k1
    real(dp), intent(in) :: u0(3),v1(3)
    real(dp), intent(out) :: u1(3)
    real(dp) :: n0,n1,u0v1(3),u0nv1(3),nu0nv1,u1v1(3),u1nv1(3)

    n0=this%indices(k0)
    n1=this%indices(k1)

    ! Decompose u0 along v1 and v1*
    u0v1=dot_product(u0,v1)
    u0nv1=u0-u0v1*v1
    nu0nv1=norm2(u0nv1)
    if (nu0nv1<=epsilon(nu0nv1)) then
       ! Orthogonal
       u1=u0
       return
    end if
    u1v1=u0v1
    u1nv1=n0/n1*u0nv1
    u1=u1v1*v1+u1nv1*u0nv1/nu0nv1
    u1=u1/norm2(u1)
  end subroutine snde_interact
    
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

  subroutine ray_trace(k0,p0,u0,bnd,inter,k_final,p_final,u_final)
    integer, intent(in) :: k0,k_final
    real(dp), intent(in) :: p0(3),u0(3)
    real(dp), intent(out) :: p_final(3),u_final(3)
    integer :: k1,k2
    real(dp) :: p1(3),u1(3),v2(3),p2(3),u2(3)
    class(boundary) :: bnd
    class(interaction) :: inter

    p1=p0
    k1=k0
    u1=u0
    
    do
       write (*,'(I2," P",3(1X,F8.4)," U",3(1X,ES11.3))') k1,p1,u1
       if (k1==k_final) then
          p_final=p1
          u_final=u1
          exit
       end if
          
       call bnd%next(k1,p1,u1,k2,p2,v2)
       if (k2<0) then
          write (*,'(A)') "STOP"
          exit
       end if
       call inter%interact(k1,u1,k2,v2,u2)
       p1=p2
       k1=k2
       u1=u2
    end do
    p_final=p1
    u_final=u1
  end subroutine ray_trace

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

  subroutine test_cylref
    type(cylindrical) :: cyl
    type(snell_descartes) :: snl
    integer, parameter :: n=5
    real(dp) :: p0(3),u0(3),p1(3),u1(3)
    integer :: k0,k1

    allocate(cyl%rs(n),snl%indices(n))
    cyl%rs=[10.0,5.0,3.0,2.0,0.5]
    snl%indices=[1.0,1.3,1.1,1.25,1.2]
    p0=[-20.0,0.5,0.123]
    u0=[0.22,0.0,0.1]
    u0=u0/norm2(u0)
    k0=1
    call ray_trace(k0,p0,u0,cyl,snl,n,p1,u1) 
  end subroutine test_cylref
end module profcov

program profcov_p
  use profcov
  implicit none

  call test_cylref
end program profcov_p
