module timescheme_mod

use params_mod,    only : params_t
use stvec_mod,     only : stvec_abstract_t
use operator_mod,  only : operator_t

implicit none

type, abstract, public :: timescheme_abstract_t
    contains
    procedure(step),  deferred :: step
end type timescheme_abstract_t

abstract interface
    subroutine step(this, f0, params, dt)
        import stvec_abstract_t!, model_parameters_abstract_t
        import timescheme_abstract_t
        import params_t

        class(timescheme_abstract_t),       intent(inout) :: this
        class(stvec_abstract_t),            intent(inout) :: f0
        type(params_t),                     intent(in)    :: params
        real(kind=8),                       intent(in)    :: dt
    end subroutine step
end interface

type, extends(timescheme_abstract_t) :: rk4_t
    class(operator_t), allocatable :: oper
    contains
    procedure step => step_rk4
end type rk4_t

type, extends(timescheme_abstract_t) :: rk4opt_t
    class(operator_t), allocatable :: oper
    class(stvec_abstract_t), allocatable :: k1,k2,k3,k4,y
    contains
    procedure step => step_rk4opt
end type rk4opt_t

contains

subroutine init_rk4opt(ts_rk4,oper,f0)
    class(timescheme_abstract_t), allocatable, intent(out) :: ts_rk4
    class(operator_t),       intent(in)  :: oper
    class(stvec_abstract_t), intent(in)  :: f0

    allocate(rk4opt_t :: ts_rk4)
    select type(ts_rk4)
    class is (rk4opt_t)
        ts_rk4%oper = oper
        allocate(ts_rk4%k1,ts_rk4%k2,ts_rk4%k3,ts_rk4%k4,ts_rk4%y,source=f0)
    end select
end subroutine init_rk4opt

subroutine step_rk4(this, f0, params, dt)
    class(rk4_t),            intent(inout) :: this
    class(stvec_abstract_t), intent(inout) :: f0
    type(params_t),          intent(in)    :: params
    real(kind=8),            intent(in)    :: dt

    class(stvec_abstract_t), allocatable :: k1,k2,k3,k4,y

    allocate(k1,k2,k3,k4,y,mold=f0)

    k1 = this%oper%act_fun(f0, params)
    y = f0+0.5_8*dt*k1
    k2 = this%oper%act_fun(y, params)
    y = f0+0.5_8*dt*k2
    k3 = this%oper%act_fun(y, params)
    y = f0+dt*k3
    k4 = this%oper%act_fun(y, params)

    f0 = f0+(k1+2._8*k2+2._8*k3+k4)*(dt/6._8)

    deallocate(k1,k2,k3,k4,y)
end subroutine step_rk4

subroutine step_rk4opt(this, f0, params, dt)
    class(rk4opt_t),         intent(inout) :: this
    class(stvec_abstract_t), intent(inout) :: f0
    type(params_t),          intent(in)    :: params
    real(kind=8),            intent(in)    :: dt

    call this%oper%act_sub(this%k1,f0, params)
    this%y = f0
    call this%y%linear_comb(this%k1,1._8,0.5_8*dt)
    call this%oper%act_sub(this%k2,this%y, params)
    this%y = f0
    call this%y%linear_comb(this%k2,1._8,0.5_8*dt)
    call this%oper%act_sub(this%k3,this%y, params)
    this%y = f0
    call this%y%linear_comb(this%k3,1._8,dt)
    call this%oper%act_sub(this%k4,this%y, params)

    call f0%linear_comb(this%k1,1._8,dt/6._8)
    call f0%linear_comb(this%k2,1._8,dt/3._8)
    call f0%linear_comb(this%k3,1._8,dt/3._8)
    call f0%linear_comb(this%k4,1._8,dt/6._8)

end subroutine step_rk4opt

end module timescheme_mod
