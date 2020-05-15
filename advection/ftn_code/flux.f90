module flux_mod
    implicit none

    type, abstract :: flux_abstract_t
        integer(kind=4) :: hw !width of halo stencil
        contains
        procedure(flux_fun_ifc), deferred :: calc_flux
    end type flux_abstract_t

    interface
        real(kind=8) function flux_fun_ifc(this, q, u) result(fl)
            import flux_abstract_t
            class(flux_abstract_t), intent(in) :: this
            real(kind=8),           intent(in) :: q(-this%hw+1:this%hw)
            real(kind=8),           intent(in) :: u
        end function flux_fun_ifc
    end interface

    type, extends(flux_abstract_t) :: up1_t
        contains
        procedure calc_flux=>calc_flux_up1
    end type up1_t
    type, extends(flux_abstract_t) :: up4_t
        contains
        procedure calc_flux=>calc_flux_up4
    end type up4_t

    contains

    type(up1_t) function init_flux_up1() result(up1_flux)
        up1_flux = up1_t(hw=1)
    end function init_flux_up1

    type(up4_t) function init_flux_up4() result(up4_flux)
        up4_flux = up4_t(hw=3)
    end function init_flux_up4


    real(kind=8) function calc_flux_up1(this, q, u) result(fl)
        class(up1_t),   intent(in) :: this
        real(kind=8), intent(in) :: q(-this%hw+1:this%hw)
        real(kind=8), intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = za1*q(0)+za2*q(1)
    end function calc_flux_up1

    real(kind=8) function calc_flux_up4(this, q, u) result(fl)
        class(up4_t),   intent(in) :: this
        real(kind=8), intent(in) :: q(-this%hw+1:this%hw)
        real(kind=8), intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = (za1*(3.0_8*q(1)+13.0_8*q(0)-5.0_8*q(-1)+q(-2))+ &
              za2*(3.0_8*q(0)+13.0_8*q(1)-5.0_8*q( 2)+q( 3)))/12.0_8
    end function calc_flux_up4

    real(kind=8) function up4f(hw, q, u) result(fl)
        integer(kind=4), intent(in) :: hw
        real(kind=8),    intent(in) :: q(-hw+1:hw)
        real(kind=8),    intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = (za1*(3.0_8*q(1)+13.0_8*q(0)-5.0_8*q(-1)+q(-2))+ &
              za2*(3.0_8*q(0)+13.0_8*q(1)-5.0_8*q( 2)+q( 3)))/12.0_8
    end function up4f

end module flux_mod
