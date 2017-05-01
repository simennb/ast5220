program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  implicit none

  ! Initialize time grids
  write(*,*) 'initialize_time_mod'
  call initialize_time_mod
  write(*,*) 'initialize_rec_mod'
  call initialize_rec_mod
  write(*,*) 'initialize_perturbation_eqns'
  call initialize_perturbation_eqns
  write(*,*) 'integrate_perturbation_eqns'
  call integrate_perturbation_eqns

  ! Output to file desired quantities here
  !write(*,*) 'David is not here. We can do anything we want :-)'
  call write_to_file_mk3

end program cmbspec
