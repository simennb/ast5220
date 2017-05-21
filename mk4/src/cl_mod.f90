module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline, n_hires, n_cl, m ! HEY CHRISTMAS TREE
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand
    real(dp),     allocatable, dimension(:,:)     :: j_l, j_l2
    real(dp),     allocatable, dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     allocatable, dimension(:)       :: k, x
    real(dp),     allocatable, dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable, dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
    real(dp),     allocatable, dimension(:)       :: x_hires, k_hires, cls_hires, ls_hires, x_lores

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Task: Get source function from evolution_mod
    n_hires = 5000
    allocate(x_hires(n_hires))
    allocate(k_hires(n_hires))
    allocate(S(n_hires,n_hires))

    call get_hires_source_function(k_hires, x_hires, S)

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))
    allocate(integrand(n_hires))

    ! Initializing z array
    do i = 1, n_spline
       z_spline(i) = (i-1.d0)*3500.d0/(n_spline-1.d0)
    end do

    ! Checking for binary file
    filename = 'j_l.bin'
    inquire(file=filename, exist=exist)
    if (exist) then
       open(10, form='unformatted', file=filename)
       read(10) j_l, j_l2
       close(10)
    else
       ! Computing spherical Bessel functions
       do l = 1, l_num
          do i = 1, n_spline
             if (z_spline(i) > 0.d0) then ! BECAUSE I'M GOING TO SAY PLEASE
                call sphbes(ls(l), z_spline(i), j_l(i,l))
             end if ! YOU HAVE NO RESPECT FOR LOGIC
          end do
       end do
       ! Splining Bessel functions
       do l = 1, l_num
          call spline(z_spline, j_l(:,l), 1.d30, 1.d30, j_l2(:,l))
       end do
       ! Write to file
       open(10, form='unformatted', file=filename)
       write(10) j_l, j_l2
       close(10)
    end if

    ! write integrand to file to check
    open(1,file='../results/integrand.dat')
    j = 2916
    l = 17
    do i = 1, n_hires
       integrand(i) = S(i,j)*splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(j)*(get_eta(0.d0)-get_eta(x_hires(i)))) 
       write(1,*) integrand(i)
    end do

    ! Overall task: Compute the C_l's for each given l
    allocate(Theta(l_num,n_hires))
    allocate(int_arg(n_hires)) ! i dont know what im doing
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(x_lores(n_hires/10)) ! creating low-res x grid for fast Theta-integration
    do l = 1, l_num

       ! Task: Compute the transfer function, Theta_l(k)
       do j = 1, n_hires
          ! Compute integrand for current k
          do i = 1, n_hires/10
             m = 1 + (i-1)*(n_hires-1)/(n_hires/10-1) ! i just wanna speed up integration, ok
             x_lores(i) = x_hires(m)
             integrand(i)=S(m,j)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(j)*(get_eta(0.d0)-get_eta(x_hires(m)))) 
          end do

          ! integrate with trapezoidal, maybe improve later
          call trapz(x_lores, integrand(1:500), Theta(l,j))

       end do

       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
       do j = 1, n_hires
          int_arg(j) = (c*k_hires(j)/H_0)**(n_s-1.d0)*Theta(l,j)**2/k_hires(j)
       end do
       call trapz(k_hires, int_arg, integral)

       ! Task: Store C_l in an array. Optionally output to file
       cls(l) = integral*ls(l)*(ls(l)+1.d0)/(2.d0*pi)

       write(*,*) ls(l), cls(l)

       ! write C_l integrand to file
       ! TODO: fix why plot is wrong
       if (ls(l) == 2) then
          call write_cl_int(ls(l), c*k_hires/H_0, Theta(l,:)**2/k_hires
       else if (ls(l) == 50) then
          call write_cl_int(ls(l), c*k_hires/H_0, Theta(l,:)**2/k_hires
       else if (ls(l) == 200) then
          call write_cl_int(ls(l), c*k_hires/H_0, Theta(l,:)**2/k_hires
       else if (ls(l) == 500) then
          call write_cl_int(ls(l), c*k_hires/H_0, Theta(l,:)**2/k_hires
       else if (ls(l) == 800) then
          call write_cl_int(ls(l), c*k_hires/H_0, Theta(l,:)**2/k_hires
       else if (ls(l) == 1200) then
          call write_cl_int(ls(l), c*k_hires/H_0, Theta(l,:)**2/k_hires
       end if

    end do

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    n_cl = ls(l_num)-1 ! 1999
    allocate(ls_dp(l_num))
    allocate(ls_hires(ls(l_num)-1)) ! 1199 different l's, could use n_hires for 5k values, but maybe overkill
    allocate(cls_hires(ls(l_num)-1))
    
    do l = 1, l_num
       ls_dp(l) = ls(l) ! spline/splint requires dp array as input
    end do
    
    call spline(ls_dp, cls, 1.d30, 1.d30, cls2) ! spline
    do i = 1, ls(l_num)-1
       ls_hires(i)  = ls_dp(1) + (i-1.d0)*(ls_dp(l_num)-ls_dp(1))/(ls_dp(l_num)-2.d0) ! -2 cause reasons
       cls_hires(i) = splint(ls_dp, cls, cls2, ls_hires(i))
    end do

    ! Write C_l's to file
    open(20,file='../results/cl.dat') ! maybe add model parameters to file name for when changing
    do i = 1, ls(l_num)-1
       write(20,'(2E20.8)') ls_hires(i), cls_hires(i)
    end do
    close(20)

  end subroutine compute_cls


  subroutine trapz(x,y,integral)
    implicit none
    real(dp), dimension(:), intent(in)  :: x,y
    real(dp),               intent(out) :: integral
    integer(i4b)                        :: n, i

    if (size(x) .ne. size(y)) then
       write(*,*) 'x and y does not have same shape'
       return
    end if

    integral = 0.d0
    n = size(x)
    do i=1, n-1
       integral = integral + (x(i+1)-x(i))*(y(i+1)+y(i))/2.d0
    end do

  end subroutine trapz


  subroutine write_cl_int(l, ckH0, cls_int)
    implicit none
    integer(i4b),            intent(in) :: l
    real(dp), dimension(:), intent(in) :: ckH0, cls_int
    integer(i4b)                       :: n, i
    character(len=128)                 :: filename, str1, str2, str3
    if (size(ckH0) .ne. size(cls_int)) then
       write(*,*) 'ckH0 and cls_int does not have same shape'
       return
    end if
    n = size(ckH0)

    ! Creating filename (cls_int_<l>.dat)
    str1 = '../results/cls_int_'
    write(str2,'(I5.5)') l
    str3 = '.dat'
    filename = trim(str1)//trim(str2)//trim(str3)
    
    open(30,file=filename)
    do i = 1, n
       write(30,'(2E20.8E4)') ckH0(i), cls_int(i)
    end do
    close(30)
  end subroutine write_cl_int
  
end module cl_mod
