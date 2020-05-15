module operator_mod

    use stvec_mod,  only: stvec_abstract_t
    use params_mod, only: params_t
    use flux_mod,   only: flux_abstract_t

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
            class(params_t),         intent(in)  :: params
        end subroutine act_sub_ifc
    end interface

    type, extends(operator_t) :: adv_oper_t
        class(flux_abstract_t), allocatable :: flux
        contains
        procedure act_fun => adv_oper_fun
        procedure act_sub => adv_oper_sub
    end type adv_oper_t

contains

    type(adv_oper_t) function init_adv_operator(namfname) result(op)
        use flux_mod, only: flux_abstract_t, init_flux_up1, init_flux_up4
        character(*)   :: namfname
        character(256) :: flux_scheme_name="up1"
        namelist /adv/ flux_scheme_name

        open(117,file=namfname, form="formatted")
        read(117,adv)
        close(117)

        if(trim(flux_scheme_name) == "up1") then
            op%flux = init_flux_up1()
        else if(trim(flux_scheme_name) == "up4") then
            op%flux = init_flux_up4()
        else
            print *, "Error: unknown flux scheme: "//flux_scheme_name
            stop
        end if
    end function init_adv_operator

    subroutine adv_oper_sub(this, fout, fin, params)

        use stvec_mod, only: stvec_t

        class(adv_oper_t),       intent(in)    :: this
        class(stvec_abstract_t), intent(inout) :: fout
        class(stvec_abstract_t), intent(in)    :: fin
        class(params_t),         intent(in)  :: params

        integer(kind=4) N

        select type(fout)
        class is (stvec_t)
        select type (fin)
        class is (stvec_t)
            N = fin%N
            call flux_conv(fout, fin, params, this%flux)
        class default
            print *, "wrong input stvec type in adv_oper_sub"
        end select
        class default
            print *, "wrong out stvec type in adv_oper_sub"
        end select

    end subroutine adv_oper_sub

    function adv_oper_fun(this, fin, params) result(fout)

        use stvec_mod, only: stvec_t

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

    subroutine flux_conv(fout, fin, params, flux)

        use stvec_mod, only: stvec_t
        use flux_mod,  only: flux_abstract_t, up4f

        class(stvec_t),         intent(inout) :: fout
        class(stvec_t),         intent(in)    :: fin
        class(params_t),        intent(in)    :: params
        class(flux_abstract_t), intent(in)    :: flux

        real(kind=8)     :: q(-flux%hw+1:params%N+flux%hw, -flux%hw+1:params%N+flux%hw)
        real(kind=8)     :: flx(0:params%N,params%N)
        real(kind=8)     :: fly(params%N,0:params%N)
        integer(kind=4)  :: N, i, j, hw

        N = params%N
        hw = flux%hw
        call periodic_bc(q,fin%p,N,hw)

        do j=1,N
            do i=0,N
                flx(i,j) = flux%calc_flux(q(i-hw+1:i+hw,j),params%u(i,j))
                !flx(i,j) = up4f(hw,q(i-hw+1:i+hw,j),params%u(i,j))
            end do
        end do
        do j=0,N
            do i=1,N
                fly(i,j) = flux%calc_flux(q(i,j-hw+1:j+hw),params%v(i,j))
                !fly(i,j) = up4f(hw,q(i,j-hw+1:j+hw),params%v(i,j))
            end do
        end do

        do j=1,N
            do i=1,N
                fout%p(i,j) =-(flx(i,j)-flx(i-1,j)+fly(i,j)-fly(i,j-1))/params%dx
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
