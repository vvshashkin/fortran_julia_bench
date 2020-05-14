module params_mod

    type params_t
        integer(kind=4)           :: N
        real(kind=8)              :: dx
        real(kind=8), allocatable :: u(:,:)
        real(kind=8), allocatable :: v(:,:)
    end type params_t

contains
    function init_params(N) result(params)
        integer(kind=4), intent(in) :: N
        type(params_t)              :: params

        params%N  = N
        params%dx = 1._8/N

        allocate(params%u(N,N))
        allocate(params%v(N,N))
        params%u(1:N,1:N) = 1._8
        params%v(1:N,1:N) = 1._8
    end function init_params

    function init_f(params) result(f)
        type(params_t), intent(in) :: params
        real(kind=8), allocatable  :: f(:,:)

        real(kind=8) r2
        real(kind=8), parameter :: x0 = 0.5_8, y0 = 0.5_8
        real(kind=8), parameter :: d2 = 0.1**2
        integer(kind=4) i, j

        allocate(f(params%N,params%N))
        do j = 1, params%N
            do i = 1, params%N
                r2 = ((i-0.5_8)*params%dx-x0)**2+((j-0.5_8)*params%dx-y0)**2
                f(i,j) = exp(-r2/d2)
            end do
        end do
    end function init_f

end module params_mod