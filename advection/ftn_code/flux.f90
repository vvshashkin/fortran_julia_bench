module flux_mod
    implicit none

    abstract interface
        real(kind=8) function clflux(hw,q,u) result(fl)
            integer(kind=4), intent(in) :: hw
            real(kind=8),    intent(in) :: q(-hw+1:hw)
            real(kind=8),    intent(in) :: u
        end function clflux
    end interface

    contains

    real(kind=8) function up1f(hw, q, u) result(fl)
        integer(kind=4), intent(in) :: hw
        real(kind=8),    intent(in) :: q(-hw+1:hw)
        real(kind=8),    intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = za1*q(0)+za2*q(1)
    end function up1f

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
