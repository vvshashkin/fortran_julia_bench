module operator_mod

    use stvec_mod,  only: stvec_abstract_t
    use params_mod, only: params_t
    use flux_mod,   only: clflux

    implicit none

    type, abstract :: operator_t
        contains
        procedure(act_fun_ifc), deferred :: act_fun
        procedure(act_sub_ifc), deferred :: act_sub
    end type operator_t

    interface
        function act_fun_ifc(this, fin, params) result(fout)
            import operator_t
            import stvec_abstract_t
            import params_t
            class(operator_t),       intent(in)  :: this
            class(stvec_abstract_t), intent(in)  :: fin
            class(params_t),         intent(in)  :: params
            class(stvec_abstract_t), allocatable :: fout
        end function act_fun_ifc

        subroutine act_sub_ifc(this, fout, fin, params)
            import operator_t
            import stvec_abstract_t
            import params_t
            class(operator_t),       intent(in)    :: this
            class(stvec_abstract_t), intent(inout) :: fout
            class(stvec_abstract_t), intent(in)    :: fin
            class(params_t),         intent(in)    :: params
        end subroutine act_sub_ifc
    end interface

    integer(kind=4), parameter :: maxhw = 4

    type, extends(operator_t) :: adv_oper_t
        integer(kind=4)                     :: hw
        procedure(clflux), pointer, nopass  :: fluxfun
        character(:), allocatable           :: flux_scheme_name
        contains
        procedure :: act_fun => adv_oper_fun
        procedure :: act_sub => adv_oper_sub
    end type adv_oper_t

contains

    type(adv_oper_t) function init_adv_operator(namfname) result(op)
        use flux_mod, only: up4f, up1f
        character(*)   :: namfname
        character(256) :: flux_scheme_name="up1"
        namelist /adv/ flux_scheme_name

        open(117,file=namfname, form="formatted")
        read(117,adv)
        close(117)

        if(trim(flux_scheme_name) == "up1") then
            op%fluxfun => up1f
            op%hw = 1
        else if(trim(flux_scheme_name) == "up4") then
            op%fluxfun => up4f
            op%hw = 3
        else
            print *, "Error: unknown flux scheme: "//flux_scheme_name
            stop
        end if

        op%flux_scheme_name = trim(flux_scheme_name)

    end function init_adv_operator

    subroutine adv_oper_sub(this, fout, fin, params)

        use stvec_mod, only: stvec_t
        use flux_mod,  only: up1f,up4f

        class(adv_oper_t),       intent(in)    :: this
        class(stvec_abstract_t), intent(inout) :: fout
        class(stvec_abstract_t), intent(in)    :: fin
        class(params_t),         intent(in)  :: params

        integer(kind=4) N

        select type(fout)
        class is (stvec_t)
        select type (fin)
        class is (stvec_t)
#ifdef ADV_OPT
            if(this%flux_scheme_name == "up1") then
                call flux_conv(fout%p, fin%p, params, up1f, 1)
            else if(this%flux_scheme_name == "up4") then
                call flux_conv(fout%p, fin%p, params, up4f, 3)
            end if
#elif HARD_CODE
            call flux_conv(fout%p, fin%p, params, this%fluxfun, 3)
#else
            call flux_conv(fout%p, fin%p, params, this%fluxfun, this%hw)
#endif
        class default
            print *, "wrong input stvec type in adv_oper_sub"
        end select
        class default
            print *, "wrong out stvec type in adv_oper_sub"
        end select

    end subroutine adv_oper_sub

    function adv_oper_fun(this, fin, params) result(fout)

        use stvec_mod, only: stvec_t
        use flux_mod,  only: up4f

        class(adv_oper_t),       intent(in)  :: this
        class(stvec_abstract_t), intent(in)  :: fin
        class(params_t),         intent(in)  :: params
        class(stvec_abstract_t), allocatable :: fout

        integer(kind=4) :: N

        allocate(stvec_t :: fout)

        select type(fout)
        class is (stvec_t)
        select type (fin)
        class is (stvec_t)
            N = fin%N
            allocate(fout%p(1:N,1:N))
            fout%N = N
            call this%act_sub(fout,fin, params)
        class default
            print *, "wrong input stvec type in adv_oper_fun"
        end select
        class default
            print *, "wrong out stvec type in adv_oper_fun"
        end select
    end function adv_oper_fun

    subroutine flux_conv(fout, fin, params, fluxf, hw)

        class(params_t),        intent(in)    :: params
        real(kind=8),           intent(inout) :: fout(params%N,params%N)
        real(kind=8),           intent(in)    :: fin(params%N,params%N)
        procedure(clflux)                     :: fluxf
        integer(kind=4),        intent(in)    :: hw

        real(kind=8)     :: q(-hw+1:params%N+hw, -hw+1:params%N+hw)
        real(kind=8)     :: flx(0:params%N,params%N)
        real(kind=8)     :: fly(params%N,0:params%N)
        integer(kind=4)  :: N, i, j
#ifdef HARD_CODE
        real(kind=8)     :: za1, za2
#endif

        N = params%N
        call periodic_bc(q,fin,N,hw)

        do j=1,N
            do i=0,N
#ifndef HARD_CODE
                flx(i,j) = fluxf(hw,q(i-hw+1:i+hw,j),params%u(i,j))
#else
                za1 = .5_8+sign(.5_8,params%u(i,j))
                za2 = 1._8-za1
                flx(i,j) = params%u(i,j)*(za1*(3.0_8*q(i+1,j)+13.0_8*q(i  ,j)-5.0_8*q(i-1,j)+q(i-2,j))+ &
                                          za2*(3.0_8*q(i  ,j)+13.0_8*q(i+1,j)-5.0_8*q(i+2,j)+q(i+3,j)))/12.0_8
#endif
            end do
        end do
        do j=0,N
            do i=1,N
#ifndef HARD_CODE
                fly(i,j) = fluxf(hw,q(i,j-hw+1:j+hw),params%v(i,j))
#else
                za1 = .5_8+sign(.5_8,params%v(i,j))
                za2 = 1._8-za1
                fly(i,j) = params%v(i,j)*(za1*(3.0_8*q(i,j+1)+13.0_8*q(i,j  )-5.0_8*q(i,j-1)+q(i,j-2))+ &
                                          za2*(3.0_8*q(i,j  )+13.0_8*q(i,j+1)-5.0_8*q(i,j+2)+q(i,j+3)))/12.0_8
#endif
            end do
        end do

        do j=1,N
            do i=1,N
                fout(i,j) =-(flx(i,j)-flx(i-1,j)+fly(i,j)-fly(i,j-1))/params%dx
            end do
        end do
    end subroutine flux_conv

    subroutine periodic_bc(q,f,N,m)
        real(kind=8),    intent(out) :: q(-m+1:N+m,-m+1:N+m)
        real(kind=8),    intent(in)  :: f(1:N,1:N)
        integer(kind=4), intent(in)  :: N, m

        q(1:N,1:N)     = f(1:N,1:N)
        q(1:N,-m+1:0)  = q(1:N,N-m+1:N)
        q(1:N,N+1:N+m) = q(1:N,1:m)
        q(-m+1:0,-m+1:N+m)  = q(N-m+1:N,-m+1:N+m)
        q(N+1:N+m,-m+1:N+m)  = q(1:m,-m+1:N+m)
    end subroutine periodic_bc

end module operator_mod
