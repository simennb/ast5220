module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),                private :: x_init   = log(a_init)
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! increasing x
  real(dp), allocatable, dimension(:)     :: x_t2

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), allocatable, dimension(:),   intent(inout) :: k, x
    real(dp), allocatable, dimension(:,:), intent(inout) :: S

    integer(i4b) :: i, j, n_hires, ii
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi, xc, kc
    real(dp), allocatable, dimension(:,:)     :: S_lores
    real(dp), allocatable, dimension(:,:,:,:) :: S_coeff

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    allocate(S_lores(1000, n_k))
    allocate(S_coeff(4,4,1000,n_k))  ! okay.....

    n_hires = size(x)
    do i = 1, n_hires  ! create hires x,k arrays
       x(i) = x_t2(501) + (i-1.d0)*(0.d0 - x_t2(501))/(n_hires-1.d0) ! change 501 to 1 if i revert to x_t
       k(i) = k_min + (k_max - k_min)*((i-1.d0)/(n_hires-1.d0))!**2
    end do
    
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    do i = 1, n_t!500!n_t  !  STICK AROUND
       ii    = i !500+i ! why did i ever expand the x-grid :/
       xc    = x_t2(i)
       g     = get_g(xc)
       dg    = get_dg(xc)
       ddg   = get_ddg(xc)
       tau   = get_tau(xc)
       dt    = get_dtau(xc)
       ddt   = get_ddtau(xc)
       H_p   = get_H_p(xc)
       dH_p  = get_dH_p(xc)
       
       do j = 1, n_k
          kc    = ks(j)
          Pi    = Theta(ii,2,j)
          dPi   = dTheta(ii,2,j)

          ! new attempt as per cmbspec_eqns document (with typos fixed) D:
          ddPi  = 2.d0*c*kc/(5.d0*H_p)*(-dH_p/H_p*Theta(ii,1,j) + dTheta(ii,1,j)) + 3.d0/10.d0*(ddt*Pi+dt*dPi) - 3.d0*c*kc/(5.d0*H_p)*(-dH_p/H_p*Theta(ii,3,j) + dTheta(ii,3,j))

          ddHH_p = H_0**2/2.d0*( (Omega_b+Omega_m)*exp(-xc) + 4.d0*Omega_r*exp(-2.d0*xc) + 4.d0*Omega_lambda*exp(2.d0*xc) )

          S_lores(i,j) = g*(Theta(ii,0,j)+Psi(ii,j)+Pi/4.d0) + exp(-tau)*(dPsi(ii,j)-dPhi(ii,j)) - 1.d0/(c*kc)*(dH_p*g*v_b(ii,j)+H_p*dg*v_b(ii,j)+H_p*g*dv_b(ii,j)) + 3.d0/(4.d0*c**2*kc**2)*(ddHH_p*g*Pi + 3.d0*H_p*dH_p*(dg*Pi+g*dPi) + H_p**2*(ddg*Pi+2.d0*dg*dPi+g*ddPi))
       end do
    end do  !  CHILL

    !   2) Then spline this function with a 2D spline
    call splie2_full_precomp(x_t2, ks, S_lores, S_coeff)

    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays
    do i=1, n_hires
       do j=1, n_hires
          S(i,j) = splin2_full_precomp(x_t2, ks, S_coeff, x(i), k(j))
       end do
    end do

  end subroutine get_hires_source_function


  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do i = 1, n_k
       ks(i) = k_min + (k_max - k_min)*((i-1.d0)/99.d0)**2
    end do

    ! increase size of x_t D:
    n_t = 1000
    allocate(x_t2(0:n_t))
    do i = 0, n_t
       if (i <= 500) then
          x_t2(i) = x_init + i*(-log(1.d0+1630.4d0) - x_init)/(500.d0+1.d0)
       else
          x_t2(i) = x_t(i-500)
       end if
    end do    

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1.d0
    delta(0,:)   = 3.d0/2.d0*Phi(0,:)
    delta_b(0,:) = 3.d0/2.d0*Phi(0,:)
       
    do i = 1, n_k
       v(0,i)       = c*ks(i)/(2*get_H_p(x_init))*Phi(0,i)
       v_b(0,i)     = c*ks(i)/(2*get_H_p(x_init))*Phi(0,i)
       Theta(0,0,i) = 0.5d0*Phi(0,i)
       Theta(0,1,i) = -c*ks(i)/(6.d0*get_H_p(x_init))*Phi(0,i)
       Theta(0,2,i) = -20.d0*c*ks(i)/(45.d0*get_H_p(x_init)*get_dtau(x_init))*Theta(0,1,i)
       do l = 3, lmax_int
          Theta(0,l,i) = -l/(2.d0*l+1)*c*ks(i)/(get_H_p(x_init)*get_dtau(x_init))*Theta(0,l-1,i)
       end do
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    x_init = log(a_init)
    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k

       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-5

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today

       do i = 1, n_t
          x1 = x_t2(i-1)
          x2 = x_t2(i)
          if (x2 <= x_tc) then
             call odeint(y_tight_coupling, x1, x2, eps, h1, hmin, tight_coupling_derivs, bsstep, output)             
             y(1:7) = y_tight_coupling(1:7)
             y(8)   = -20.d0*c*k_current/(45.d0*get_H_p(x2)*get_dtau(x2))*y(7)
             do l = 3, lmax_int
                y(6+l) = -l/(2.d0*l+1.d0)*c*k_current/(get_H_p(x2)*get_dtau(x2))*y(6+l-1)
             end do
          else
             ! Task: Integrate equations from tight coupling to today
             call odeint(y, x1, x2, eps, h1, hmin, derivs, bsstep, output)
          end if

          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)
          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l) 
          end do
          Psi(i,k)     = -y(5)-12.d0*H_0**2/(c**2*k_current**2*exp(2.d0*x2))*Omega_r*y(8)

          ! Task: Store derivatives that are required for C_l estimation
          !  fixed in milestone 4: actually store correct derivs in tight-coupling regime
          if (x2 <= x_tc) then
             call tight_coupling_derivs(x2, y, dydx)
             dv_b(i,k)     = dydx(4)
             dPhi(i,k)     = dydx(5)
             dTheta(i,0,k) = dydx(6)
             dTheta(i,1,k) = dydx(7)
             dTheta(i,2,k) = 20.d0*c*k_current*Theta(i,1,k)/(45.d0*get_H_p(x2)*get_dtau(x2))*(get_dH_p(x2)/get_H_p(x2)+get_ddtau(x2)/get_dtau(x2)-dTheta(i,1,k)/Theta(i,1,k))
             do l = 3, lmax_int
                dTheta(i,l,k) = l*c*k_current*Theta(i,l-1,k)/((2.d0*l+1.d0)*get_H_p(x2)*get_dtau(x2))*(get_dH_p(x2)/get_H_p(x2)+get_ddtau(x2)/get_dtau(x2)-dTheta(i,l-1,k)/Theta(i,l-1,k))
             end do
          else
             call derivs(x2, y, dydx)
             dv_b(i,k)     = dydx(4)
             dPhi(i,k)     = dydx(5)
             dTheta(i,:,k) = dydx(6:6+lmax_int)
          end if
          dPsi(i,k)     = -dydx(5) - 12.d0*H_0**2/(c**2*k_current**2*exp(2.d0*x2))*Omega_r*(dydx(8)-2.d0*y(8))
       end do
    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns


  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    integer(i4b)          :: i, n_test
    real(dp)              :: get_tight_coupling_time, x, x_start_rec

    n_test = 1000
    x_start_rec = -log(1.d0+1630.4d0)
    x = x_init
    do i = 1, n_test
       x = x + (x_start_rec-x_init)/(n_test-1)
       if ((abs(c*k/(get_H_p(x)*get_dtau(x))) >= 0.1d0 .and. abs(get_dtau(x))<=10.d0) .or. ( x >= x_start_rec)) then
          ! (abs(get_dtau(x)) <= 10.d0 .and.
          get_tight_coupling_time = x
          exit
       end if
    end do

  end function get_tight_coupling_time


  ! Subroutines for calculating the derivatives for integration
  subroutine tight_coupling_derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: Psi, Theta_2, H_p, a, R, q, dtau, ckH_p, dH_p, ddtau

    ! Initializing some variables to be using in computing derivatives
    H_p     = get_H_p(x)
    dH_p    = get_dH_p(x)
    a       = exp(x)
    dtau    = get_dtau(x)
    ddtau   = get_ddtau(x)
    R       = 4.d0*Omega_r/(3.d0*Omega_b*a)
    ckH_p   = c*k_current/H_p
    Theta_2 = -20.d0*ckH_p/(45.d0*dtau)*y(7)
    Psi     = -y(5) - 12.d0*(H_0/(c*k_current*a))**2*Omega_r*Theta_2

    ! y(5)  =  Phi(0,k)
    dydx(5) = Psi - ckH_p**2/3.d0*y(5) + H_0**2/(2.d0*H_p**2)*( Omega_m/a*y(1) + Omega_b/a*y(2) + 4.d0*Omega_r/a**2*y(6) )

    ! y(1)  =  delta(0,k)
    dydx(1) = ckH_p*y(3) - 3.d0*dydx(5)

    ! y(2)  =  delta_b(0,k)
    dydx(2) = ckH_p*y(4) - 3.d0*dydx(5)

    ! y(3)  =  v(0,k)
    dydx(3) = -y(3) - ckH_p*Psi

    ! y(6)  =  Theta(0,0,k)
    dydx(6) = -ckH_p*y(7) - dydx(5)

    ! y(4)  =  v_b(0,k)
    q = ( -( (1.d0-2.d0*R)*dtau + (1.d0+R)*ddtau )*(3.d0*y(7)+y(4)) - ckH_p*Psi + (1.d0-dH_p/H_p)*ckH_p*(-y(6)+2.d0*Theta_2 ) - ckH_p*dydx(6) )/( (1.d0+R)*dtau+dH_p/H_p - 1.d0 )
    dydx(4) = 1.d0/(1.d0+R)*(-y(4)-ckH_p*Psi+R*(q+ckH_p*(-y(6)+2.d0*Theta_2)-ckH_p*Psi))
    
    ! y(7)  =  Theta(0,1,k)
    dydx(7) = 1.d0/3.d0*(q-dydx(4))

  end subroutine tight_coupling_derivs

  subroutine derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: Psi, H_p, a, R, dtau, ckH_p, eta
    integer(i4b)                        :: l

    ! Initializing some variables to be using in computing derivatives
    H_p     = get_H_p(x)
    eta     = get_eta(x)
    a       = exp(x)
    dtau    = get_dtau(x)
    ckH_p   = c*k_current/H_p
    R       = 4.d0*Omega_r/(3.d0*Omega_b*a)
    Psi     = -y(5) - 12.d0*H_0**2/(c**2*k_current**2*a**2)*Omega_r*y(8)

    ! y(5)  =  Phi(0,k)
    dydx(5) = Psi - ckH_p**2/3.d0*y(5) + H_0**2/(2.d0*H_p**2)*( Omega_m/a*y(1) + Omega_b/a*y(2) + 4.d0*Omega_r/a**2*y(6) )

    ! y(1)  =  delta(0,k)
    dydx(1) = ckH_p*y(3) - 3.d0*dydx(5)

    ! y(2)  =  delta_b(0,k)
    dydx(2) = ckH_p*y(4) - 3.d0*dydx(5)

    ! y(3)  =  v(0,k)
    dydx(3) = -y(3) - ckH_p*Psi

    ! y(4)  =  v_b(0,k)
    dydx(4) = -y(4) - ckH_p*Psi + dtau*R*(3.d0*y(7)+y(4))

    ! y(6)  =  Theta(0,0,k)
    dydx(6) = -ckH_p*y(7) - dydx(5)
    
    ! y(7)  =  Theta(0,1,k)
    dydx(7) = ckH_p/3.d0*y(6) - 2.d0*ckH_p/3.d0*y(8) + ckH_p/3.d0*Psi + dtau*( y(7) + 1.d0/3.d0*y(4) )

    ! y(l)  =  Theta(0,l,k) , 2 <= l < lmax_int
    do l = 2, lmax_int-1
       dydx(6+l) = l*ckH_p/(2.d0*l+1.d0)*y(6+l-1) - (l+1.d0)*ckH_p/(2.d0*l+1.d0)*y(6+l+1) + dtau*(y(6+l) - 1.d0/10.d0*y(6+l)*abs(l==2))
    end do

   ! dydx(8) = 2.d0*ckH_p/5.d0*y(7) - 3.d0*ckH_p/5.d0*y(9) + dtau*(y(8) - 1.d0/10.d0*y(8))
    !do l = 3, lmax_int-1
     !  dydx(6+l) = l*ckH_p/(2.d0*l+1.d0)*y(6+l-1) - (l+1.d0)*ckH_p/(2.d0*l+1.d0)*y(6+l+1) + dtau*(y(6+l))
    !end do
    
    ! y(l)  =  Theta(0,l,k) , l = lmax_int
    ! using that fortran for some reason increments l once more as do loop ends
    dydx(6+l) = ckH_p*y(6+l-1) - c*(l+1.d0)/(H_p*eta)*y(6+l) + dtau*y(6+l)

  end subroutine derivs

  ! ----------- Milestone 3 - Write to file -----------
  subroutine write_to_file_mk3
    use healpix_types
    implicit none
    integer(i4b)                             :: i
    integer(i4b), allocatable, dimension(:)  :: kw
    allocate(kw(n_k))
    ! Need to pick 6 values of k to show the three regimes
    kw(1:6) = (/ 1, 12, 30, 40, 85, 100 /)
!    kw(1:6) = (/ 1, 10, 30, 50, 75, 100 /)

    open(1, file='../results/phi.dat')
    open(2, file='../results/psi.dat')
    open(3, file='../results/delta.dat')
    open(4, file='../results/delta_b.dat')
    open(5, file='../results/v.dat')
    open(6, file='../results/v_b.dat')
    open(7, file='../results/theta.dat') ! l=0
    open(8, file='../results/k_values.dat')

    do i=1, 6
       write(8,*) ks(i), kw(i)
    end do

    do i=1, n_t
       write(1,'(7F20.7)') x_t2(i),Phi(i,kw(1)),Phi(i,kw(2)),Phi(i,kw(3)),Phi(i,kw(4)),Phi(i,kw(5)),Phi(i,kw(6))
       write(2,'(7F20.7)') x_t2(i),Psi(i,kw(1)),Psi(i,kw(2)),Psi(i,kw(3)),Psi(i,kw(4)),Psi(i,kw(5)),Psi(i,kw(6))
       write(3,'(7F20.7)') x_t2(i),delta(i,kw(1)),delta(i,kw(2)),delta(i,kw(3)),delta(i,kw(4)),delta(i,kw(5)),delta(i,kw(6))
       write(4,'(7F20.7)') x_t2(i),delta_b(i,kw(1)),delta_b(i,kw(2)),delta_b(i,kw(3)),delta_b(i,kw(4)),delta_b(i,kw(5)),delta_b(i,kw(6))
       write(5,'(7F20.7)') x_t2(i),v(i,kw(1)),v(i,kw(2)),v(i,kw(3)),v(i,kw(4)),v(i,kw(5)),v(i,kw(6))
       write(6,'(7F20.7)') x_t2(i),v_b(i,kw(1)),v_b(i,kw(2)),v_b(i,kw(3)),v_b(i,kw(4)),v_b(i,kw(5)),v_b(i,kw(6))
       write(7,'(7F20.7)') x_t2(i),Theta(i,0,kw(1)),Theta(i,0,kw(2)),Theta(i,0,kw(3)),Theta(i,0,kw(4)),Theta(i,0,kw(5)),Theta(i,0,kw(6))
    end do

    ! close files
    do i=1, 8
       close(i)
    end do
    deallocate(kw)

  end subroutine write_to_file_mk3

end module evolution_mod
