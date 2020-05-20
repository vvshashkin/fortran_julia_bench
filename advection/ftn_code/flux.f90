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
        fl = u*(za1*q(0)+za2*q(1))
    end function up1f

    real(kind=8) function up4f(hw, q, u) result(fl)
        integer(kind=4), intent(in) :: hw
        real(kind=8),    intent(in) :: q(-hw+1:hw)
        real(kind=8),    intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = u*(za1*(3.0_8*q(1)+13.0_8*q(0)-5.0_8*q(-1)+q(-2))+ &
                za2*(3.0_8*q(0)+13.0_8*q(1)-5.0_8*q( 2)+q( 3)))/12.0_8
    end function up4f

    real(kind=8) function weno5f(hw,q,u) result(fl)
        integer(kind=4), intent(in) :: hw
        real(kind=8),    intent(in) :: q(-hw+1:hw)
        real(kind=8),    intent(in) :: u

        real(kind=8), parameter :: eps = 1e-6_8
        real(kind=8) :: za1, za2, a, b, c, d, e
        real(kind=8) :: q1, q2, q3, s1, s2, s3, w1, w2, w3

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        a = za1*q(hw-2)+za2*q(hw+3)
        b = za1*q(hw-1)+za2*q(hw+2)
        c = za1*q(hw  )+za2*q(hw+1)
        d = za1*q(hw+1)+za2*q(hw)
        e = za1*q(hw+2)+za2*q(hw-1)

        q1 = (2._8*a-7._8*b+11._8*c)/6._8
        q2 = (-b+5._8*c+2._8*d)/6._8
        q3 = (2._8*c+5._8*d-e)/6._8
        s1 = (13._8/12._8)*(a-2._8*b+c)**2 + 0.25*(a-4._8*b+3._8*c)**2
        s2 = (13._8/12._8)*(b-2._8*c+d)**2 + 0.25*(d-b)**2
        s3 = (13._8/12._8)*(c-2._8*d+e)**2 + 0.25*(3._8*c-4._8*d+e)**2

        w1 = 0.1_8/(eps+s1)**2
        w2 = 0.6_8/(eps+s2)**2
        w3 = 0.3_8/(eps+s3)**2
        fl = u*(w1*q1+w2*q2+w3*q3)/(w1+w2+w3)
    end function weno5f

end module flux_mod
