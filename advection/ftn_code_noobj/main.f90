program advection_test
    use timescheme_mod,  only : step_rk4
    use params_mod, only : init_params, init_f

    implicit none

    real(kind=8)    :: dt = 0.0005_8
    integer(kind=4) :: N = 1000
    integer(kind=4) :: Nstep = 100
    integer(kind=4) :: it
    real(kind=8)    :: dx
    real(kind=8), allocatable :: q(:,:)

    namelist /prm/ N, dt, Nstep

    open(117,file="namelist", form="formatted")
    read(117,prm)
    close(117)


    call init_params(N)
    q = init_f()
    print *, "q t=0 max/min:", maxval(q), minval(q)

    do it = 1,Nstep
        call step_rk4(q,dt)
    end do

    print *, "q(t= ", Nstep*dt, ") max/min:", maxval(q), minval(q)

end program advection_test
