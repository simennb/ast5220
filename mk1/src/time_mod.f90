module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init, x_step, crit_param

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))

    ! From z=1630 to z=614
    x_step =  (x_end_rec - x_start_rec) / (n1-1)
    do i = 1, n1
       x_t(i) = x_start_rec+(i-1)*x_step
       a_t(i) = exp(x_t(i))
    end do

    ! From z=614 to z=0
    x_step = (x_0 - x_end_rec) / (n2)
    do i = 0, n2
       x_t(n1+i) = x_end_rec+(i)*x_step
       a_t(n1+i) = exp(x_t(n1+i))
    end do

    ! Task: 1) Compute the conformal time at each eta time step    
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))
    
    x_step = (x_eta2 - x_eta1) / (n_eta-1)
    do i = 1, n_eta
       x_eta(i) = x_eta1+ (i-1)*x_step
    end do
    
!    eta(1) = a_init/(H_0*sqrt(Omega_r)) ! initial conformal time, check NO SHOULD BE BASED ON x, NOT a
    call odeint(eta(1:1), -100.d0, x_eta(1), 1.d-10, 1.d0, 0.d0, eta_derivs, bsstep, output)

    ! WHY THE FUCK IS THIS FUCKING CODE SO FUCKING CONVOLUTED AND WHY THE FUCK IS ODEINT STUPID
    do i = 2, n_eta
       eta(i) = eta(i-1)
       call odeint(eta(i:i), x_eta(i-1), x_eta(i), 1.d-10, 1.d0, 0.d0, eta_derivs, bsstep, output)
    end do
  
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    call spline(x_eta, eta, 1.d30, 1.d30, eta2)
    

    ! writing results to file
    open(1, file='../results/x_eta.dat')
    open(2, file='../results/x_H.dat')
    open(3, file='../results/x_Omegas.dat')
    do i = 1, n_eta
       write(1,*) x_eta(i), get_eta(x_eta(i))
       write(2,*) x_eta(i), get_H(x_eta(i))

       crit_param = H_0**2/get_H(x_eta(i))**2
       write(3,'(4F10.7)') Omega_b*(exp(-3*x_eta(i)))*crit_param, Omega_m*(exp(-3*x_eta(i)))*crit_param, Omega_r*(exp(-4*x_eta(i)))*crit_param, Omega_lambda*crit_param
    end do

    do i=1, 3
       close(i)
    end do

  end subroutine initialize_time_mod
  
  subroutine eta_derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    
    dydx = c/get_H_p(x)
  end subroutine

  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output

  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H, a

    a = exp(x)
    get_H = H_0*sqrt( (Omega_b+Omega_m)*a**(-3) + Omega_r*a**(-4) + Omega_lambda)

  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p

    get_H_p = exp(x)*get_H(x)

  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p, h

    ! CHECK THIS WITH HKE OR SOMEONE ELSE IF THERE IS A BETTER WAY OF DOING IT
    ! this is just temporary so that the function returns something
    h = 0.001
    get_dH_p = (get_H_p(x+h)-get_H_p(x-h))/(2.d0*h)

  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

    get_eta = splint(x_eta, eta, eta2, x_in)

  end function get_eta

end module time_mod
