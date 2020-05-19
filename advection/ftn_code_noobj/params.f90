module params_mod

implicit none

integer(kind=4) :: N
real(kind=8)    :: dx
real(kind=8), allocatable :: u(:,:), v(:,:)

contains

    subroutine init_params(Nn)
        integer(kind=4), intent(in) :: Nn

        N = Nn
        dx = 1.0_8/N
        allocate(u(0:N,1:N),v(0:N,1:N))
        u(:,:) = 1._8
        v(:,:) = 1._8

    end subroutine init_params

    function init_f() result(f)
        real(kind=8)                :: f(N,N)

        real(kind=8) r2
        real(kind=8), parameter :: x0 = 0.5_8, y0 = 0.5_8
        real(kind=8), parameter :: d2 = 0.1**2
        integer(kind=4) i, j

        do j = 1, N
            do i = 1, N
                r2 = ((i-0.5_8)*dx-x0)**2+((j-0.5_8)*dx-y0)**2
                f(i,j) = exp(-r2/d2)
            end do
        end do

    end function init_f

end module params_mod
