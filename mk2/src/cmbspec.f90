program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  implicit none

  ! Initialize time grids
  call initialize_time_mod
  call initialize_rec_mod

  ! Output to file desired quantities here
  !write(*,*) 'David is not here. We can do anything we want :-)'

end program cmbspec
