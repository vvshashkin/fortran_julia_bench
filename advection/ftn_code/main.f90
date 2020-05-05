program advection_test
    use stvec_mod,       only: stvec_t, init_stvec
    use timescheme_mod, only : timescheme_abstract_t, rk4_t, init_rk4opt
    use operator_mod,    only : operator_t, adv_oper_t

    implicit none

    type(stvec_t)                     :: f1
    class(timescheme_abstract_t), allocatable :: ts
    class(operator_t),    allocatable :: oper

    real(kind=8), parameter    :: dt = 1.5_8
    integer(kind=4), parameter :: N = 1000
    real(kind=8)               :: p(N,N)
    integer(kind=4)            :: it

    call random_number(p)
    p(1,1) = 1._8
    call init_stvec(f1,N,p)

    oper=adv_oper_t()
    !ts = rk4_t(oper)
    call init_rk4opt(ts, oper, f1)

    do it = 1,100
        call ts%step(f1,dt)
    end do

    print *, maxval(f1%p) -exp(-0.1_8*it*dt)

end program advection_test
