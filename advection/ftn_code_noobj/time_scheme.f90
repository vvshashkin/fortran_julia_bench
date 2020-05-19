module timescheme_mod

use params_mod,    only : N
use operator_mod,  only : flux_conv

implicit none

contains

    subroutine step_rk4(f0, dt)
        real(kind=8), intent(inout) :: f0(N,N)
        real(kind=8), intent(in)    :: dt

        real(kind=8) y(N,N), k1(N,N), k2(N,N), k3(N,N), k4(N,N)

        call flux_conv(k1, f0)
        y = f0+0.5_8*dt*k1
        call flux_conv(k2, y)
        y = f0+0.5_8*dt*k2
        call flux_conv(k3, y)
        y = f0+dt*k3
        call flux_conv(k4, y)

        f0 = f0+(k1+2._8*k2+2._8*k3+k4)*(dt/6._8)

    end subroutine step_rk4

end module timescheme_mod
