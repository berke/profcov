module bodhaine
  use precision
  implicit none
  integer, parameter, private :: wp=dp
contains
  ! nu      [cm^-1]
  ! p       [atm]
  ! T       [K]
  ! rho_H2O [g/m^3]
  ! X_CO2   [ppmv]
  elemental function n_air_minus_1(nu,p,T,rho_H2O,X_CO2) result(dn)
    real(wp), intent(in) :: nu,p,T,rho_H2O,X_CO2
    real(wp), parameter :: &
         C_CO2=5.40335e-7_wp, &
         A_H2O=1.9809e-10_wp, &
         B_H2O=3.1759e-21_wp, &
         d0=0.0232226_wp, &
         d1=7.147815e8_wp, &
         d2=1.32274e10_wp, &
         d3=5.029045e6_wp, &
         d4=3.932957e9_wp
    real(wp) :: nu2,dn

    nu2=nu*nu
    dn=p*(d0+d1/(d2-nu2)+d3/(d4-nu2))*(1+C_CO2*X_CO2)/T &
         - (A_H2O-B_H2O*nu2)*rho_H2O*T
  end function n_air_minus_1
end module bodhaine
