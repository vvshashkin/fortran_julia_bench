module operator_mod

    implicit none

contains

    subroutine flux_conv(fout, fin)

        use params_mod, only : N, dx, u, v

        real(kind=8),           intent(inout) :: fout(N,N)
        real(kind=8),           intent(in)    :: fin(N,N)

        integer(kind=4), parameter :: hw = 3

        real(kind=8)     :: q(-hw+1:N+hw, -hw+1:N+hw)
        real(kind=8)     :: flx(0:N,N)
        real(kind=8)     :: fly(N,0:N)
#ifdef HARD_CODE
        real(kind=8)     :: za1, za2
#endif
        integer(kind=4)  :: i, j

        call periodic_bc(q,fin,N,hw)

        do j=1,N
            do i=0,N
#ifndef HARD_CODE
                flx(i,j) = up4f(hw,q(i-hw+1:i+hw,j),u(i,j))
#else
                za1 = .5_8+sign(.5_8,u(i,j))
                za2 = 1._8-za1
                flx(i,j) = u(i,j)*(za1*(3.0_8*q(i+1,j)+13.0_8*q(i  ,j)-5.0_8*q(i-1,j)+q(i-2,j))+ &
                                   za2*(3.0_8*q(i  ,j)+13.0_8*q(i+1,j)-5.0_8*q(i+2,j)+q(i+3,j)))/12.0_8
#endif
            end do
        end do
        do j=0,N
            do i=1,N
#ifndef HARD_CODE
                fly(i,j) = up4f(hw,q(i,j-hw+1:j+hw),v(i,j))
#else
                za1 = .5_8+sign(.5_8,v(i,j))
                za2 = 1._8-za1
                fly(i,j) = v(i,j)*(za1*(3.0_8*q(i,j+1)+13.0_8*q(i,j  )-5.0_8*q(i,j-1)+q(i,j-2))+ &
                                   za2*(3.0_8*q(i,j  )+13.0_8*q(i,j+1)-5.0_8*q(i,j+2)+q(i,j+3)))/12.0_8
#endif
            end do
        end do

        do j=1,N
            do i=1,N
                fout(i,j) =-(flx(i,j)-flx(i-1,j)+fly(i,j)-fly(i,j-1))/dx
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

    pure real(kind=8) function up4f(hw, q, u) result(fl)
        integer(kind=4), intent(in) :: hw
        real(kind=8),    intent(in) :: q(-hw+1:hw)
        real(kind=8),    intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = u*(za1*(3.0_8*q(1)+13.0_8*q(0)-5.0_8*q(-1)+q(-2))+ &
                za2*(3.0_8*q(0)+13.0_8*q(1)-5.0_8*q( 2)+q( 3)))/12.0_8
    end function up4f
end module operator_mod
