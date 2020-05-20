module params_mod

    type params_t
        integer(kind=4)           :: N
        real(kind=8)              :: dx
        real(kind=8), allocatable :: u(:,:)
        real(kind=8), allocatable :: v(:,:)
        integer(kind=4)           :: wf
        integer(kind=4)           :: fdist
    end type params_t

contains
    function init_params(N) result(params)
        integer(kind=4), intent(in) :: N
        type(params_t)              :: params

        integer(kind=4) :: wf = 1, fdist = 1
        real(kind=8)    :: zw
        namelist /testcase/ wf, fdist

        open(117,file="namelist", form="formatted")
        read(117,testcase)
        close(117)

        params%N  = N
        params%dx = 1._8/N
        params%wf = wf
        params%fdist = fdist

        allocate(params%u(0:N,1:N))
        allocate(params%v(1:N,0:N))
        if(wf == 1) then
            params%u(0:N,1:N) = 1._8
            params%v(1:N,0:N) = 1._8
        else if(wf == -1) then
            params%u(0:N,1:N) = -1._8
            params%v(1:N,0:N) =  1._8
        else if(wf == 2) then
            do j=1,N
                do i=0,N
                    zw = max(0._8, &
                             1.0_8-2.0_8*sqrt((i*params%dx-0.5_8)**2+ &
                                             ((j-0.5_8)*params%dx-0.5_8)**2))
                    params%u(i,j) = ((j-0.5)*params%dx-0.5_8)*zw
                end do
            end do
            do j=0,N
                do i=1,N
                    zw = max(0._8, &
                             1.0_8-2.0_8*sqrt(((i-0.5_8)*params%dx-0.5_8)**2+ &
                                             (j*params%dx-0.5_8)**2))
                    params%v(i,j) =-((i-0.5)*params%dx-0.5_8)*zw
                end do
            end do
        else
            print *, "Error: unknown wind-field"
        end if
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
