module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function

contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop, xstep
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstop

    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    ! Task: Fill in x (rec) grid
    xstep = (xstop - xstart) / (n-1)
    do i=1, n
       x_rec(i) = xstart + (i-1)*xstep
    end do

    ! Task: Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1,n
       if (use_saha) then
          ! Use the Saha equation
          n_b = Omega_b*rho_c / ( m_H*exp(3.d0*x_rec(i)) )
          T_b = T_0 / exp(x_rec(i))

          y = 1.d0/n_b*( m_e*k_b*T_b / (2.d0*pi*hbar**2))**(3.d0/2.d0)*exp(-epsilon_0/(k_b*T_b))
          
          X_e(i) = (-y + sqrt(y**2+4*y))/2.d0

          n_e(i) = n_b*X_e(i)

          if (X_e(i) < saha_limit) use_saha = .false.
       else
          ! Use the Peebles equation
          n_b = Omega_b*rho_c / ( m_H*exp(3.d0*x_rec(i)) )

          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), 1.d-10, xstep, 0.d0, X_e_derivs, bsstep, output)
          n_e(i) = n_b*X_e(i)
       end if
    end do
    
    ! Task: Compute splined (log of) electron density function
    call spline(x_rec, log(n_e), 1.d30, 1.d30, n_e2) 

    ! Task: Compute optical depth at all grid points
    tau(n) = 0.d0
    do i = n-1, 1, -1
          tau(i) = tau(i+1)
          call odeint(tau(i:i), x_rec(i+1), x_rec(i), 1.d-10, xstep, 0.d0, tau_derivs, bsstep, output)
    end do

    ! Task: Compute splined (log of) optical depth
    call spline(x_rec, tau, 1.d30, 1.d30, tau2)

    ! Task: Compute splined second derivative of (log of) optical depth
    call spline(x_rec, tau2, 1.d30, 1.d30, tau22)

    ! Task: Compute splined visibility function
    do i = 1, n
       g(i) = -get_dtau(x_rec(i))*exp(-get_tau(x_rec(i)))
    end do
    
    call spline(x_rec, g, 1.d30, 1.d30, g2)

    ! Task: Compute splined second derivative of visibility function
    call spline(x_rec, g2, 1.d30, 1.d30, g22)


    ! ------------- Writing to file ----------------
!    open(1, file='../results/x_X_e.dat')
!    open(2, file='../results/x_tau.dat')
!    open(3, file='../results/x_g.dat')
!    
!    do i=1, n
!       write(1,*) x_rec(i), X_e(i)
!       write(2,'(4F20.7)') x_rec(i), get_tau(x_rec(i)), get_dtau(x_rec(i)), get_ddtau(x_rec(i))
!       write(3,'(4F20.7)') x_rec(i), get_g(x_rec(i)), get_dg(x_rec(i)), get_ddg(x_rec(i))
!    end do
!
!    do i=1, 3
!       close(i)
!    end do

  end subroutine initialize_rec_mod

  ! Subroutines for using odeint for Peebles' equation
  subroutine X_e_derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: n_b,T_b, phi_2, alpha_2, beta, beta_2, n_1s, lambda_a, lambda_21, C_r

    n_b = Omega_b*rho_c / ( m_H*exp(3.d0*x) )
    T_b = T_0 / exp(x)

    ! Peebles' equation with units
    phi_2     = 0.448d0*log(epsilon_0/(k_B*T_b))
    alpha_2   = 64.d0*pi/sqrt(27.d0*pi) * (hbar*alpha)**2.d0/(m_e**2.d0*c) * sqrt(epsilon_0/(k_B*T_b)) * phi_2
    beta      = alpha_2 * (m_e*k_B*T_b / (2.d0*pi*hbar**2.d0) )**(3.d0/2.d0) * exp(-epsilon_0/(k_B*T_b))

    ! in order to avoid beta_2 going to NaN when T_b gets sufficiently small
    if (T_b > 169.d0) then
       beta_2 = beta * exp(3.d0*epsilon_0/(4.d0*k_B*T_b))
    else
       beta_2 = 0.d0
    endif
       
    n_1s      = (1.d0 - y(1))*n_b
    lambda_a  = get_H(x) * (3.d0*epsilon_0/ (c*hbar))**3.d0 / ((8.d0*pi)**2.d0 * n_1s)
    lambda_21 = 8.227d0  ! can be moved outside, constant
    C_r       = (lambda_21 + lambda_a) / (lambda_21 + lambda_a + beta_2)

    dydx = C_r/get_H(x) * ( beta*(1.d0-y(1)) - n_b*alpha_2*y(1)**2.d0 )

  end subroutine X_e_derivs

  ! Tau deriv subroutine
  subroutine tau_derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    dydx = -get_n_e(x)*sigma_T*exp(x)*c/get_H_p(x)

  end subroutine tau_derivs

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    get_n_e = exp(splint(x_rec, log(n_e), n_e2, x))

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau
    get_tau = splint(x_rec, tau, tau2, x)
 !   get_tau = exp(splint(x_rec, log(tau), tau2, x))

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau
    get_dtau = splint_deriv(x_rec, tau, tau2, x)
!    get_dtau = exp(splint_deriv(x_rec, log(tau), tau2, x))

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    get_ddtau = splint(x_rec, tau2, tau22, x)
!    get_ddtau = exp(splint(x_rec, tau2, tau22, x))

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitrary x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g
    get_g = splint(x_rec, g, g2, x)

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitrary x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg
    get_dg = splint_deriv(x_rec, g, g2, x)

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitrary x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
    get_ddg = splint(x_rec, g2, g22, x)

  end function get_ddg

end module rec_mod
