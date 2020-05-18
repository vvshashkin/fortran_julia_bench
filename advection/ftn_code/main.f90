program advection_test
    use params_mod,      only : params_t, init_params, init_f
    use stvec_mod,       only : stvec_t, init_stvec
    use timescheme_mod,  only : timescheme_abstract_t, rk4_t, init_rk4opt, init_rk4
    use operator_mod,    only : operator_t, adv_oper_t, init_adv_operator

    implicit none

    type(params_t)                            :: params
    type(stvec_t)                             :: f1
    class(timescheme_abstract_t), allocatable :: ts
    class(operator_t),            allocatable :: oper

    real(kind=8)    :: dt = 0.0005_8
    integer(kind=4) :: N = 1000
    integer(kind=4) :: Nstep = 100
    character(256)  :: time_scheme = "rk4_opt"
    integer(kind=4) :: it

    namelist /prm/ N, dt, Nstep, time_scheme

    open(117,file="namelist", form="formatted")
    read(117,prm)
    close(117)

    params = init_params(N)

    call init_stvec(f1,params%N,init_f(params))
    print *, "q t=0 max/min:", maxval(f1%p), minval(f1%p)

    oper = init_adv_operator("namelist")!adv_oper_t()
    if(trim(time_scheme) == "rk4_opt") then
        call init_rk4opt(ts, oper, f1)
    else if(trim(time_scheme) == "rk4") then
        !ts = rk4_t(oper)
        call init_rk4(ts, oper, f1)
    else
        print *, "unknown time-scheme: ", trim(time_scheme)
        stop
    end if

    do it = 1,Nstep
        call ts%step(f1,params,dt)
    end do

    print *, "q(t= ", Nstep*dt, ") max/min:", maxval(f1%p), minval(f1%p) 

end program advection_test
