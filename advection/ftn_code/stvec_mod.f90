module stvec_mod

    implicit none

    type, abstract :: stvec_abstract_t
        contains
            generic                                                :: operator(*) => scalar_mult_stvec, &
                                                                                     stvec_mult_scalar
            procedure(stvec_mult_scalar_abs), deferred             :: stvec_mult_scalar
            procedure(scalar_mult_stvec_abs), deferred, pass(this) :: scalar_mult_stvec
            !sum of two state vectors
            generic                                                :: operator(+) => stvec_add_stvec
            procedure(stvec_add_stvec_abs), deferred               :: stvec_add_stvec
            !memory safe (?) assign
            generic                                                :: assignment(=) => stvec_assign_stvec
            !procedure(stvec_assign_stvec_abs), deferred, pass(lhs) :: stvec_assign_stvec
            !this line to avoid linking issues with stvec_assign_stvec_abs when
            !stvec_abstract_t is a component of some type
            procedure, pass(lhs) :: stvec_assign_stvec => stvec_assign_stvec_base

            procedure(linear_comb_ifc), deferred :: linear_comb
            procedure(smult_ifc),       deferred :: smult
    end type stvec_abstract_t

    interface
        function scalar_mult_stvec_abs(a, this) result(mult)
            import stvec_abstract_t
            class(stvec_abstract_t), intent(in)  :: this
            real(kind=8),            intent(in)  :: a
            class(stvec_abstract_t), allocatable :: mult
        end function scalar_mult_stvec_abs

        function stvec_mult_scalar_abs(this,a) result(mult)
            import stvec_abstract_t
            class(stvec_abstract_t), intent(in)  :: this
            real(kind=8),            intent(in)  :: a
            class(stvec_abstract_t), allocatable :: mult
        end function stvec_mult_scalar_abs

        function stvec_add_stvec_abs(this, other) result(summed)
            import stvec_abstract_t
            class(stvec_abstract_t), intent(in)  :: this
            class(stvec_abstract_t), intent(in)  :: other
            class(stvec_abstract_t), allocatable :: summed
        end function stvec_add_stvec_abs

        subroutine stvec_assign_stvec_abs(lhs, rhs)
            import stvec_abstract_t
            class(stvec_abstract_t), intent(inout) :: lhs
            class(stvec_abstract_t), intent(in)    :: rhs
        end subroutine stvec_assign_stvec_abs

        subroutine linear_comb_ifc(this,other,a,b)
            import stvec_abstract_t
            class(stvec_abstract_t), intent(inout) :: this
            class(stvec_abstract_t), intent(in)    :: other
            real(kind=8),            intent(in)    :: a, b
        end subroutine linear_comb_ifc
        subroutine smult_ifc(this,a)
            import stvec_abstract_t
            class(stvec_abstract_t), intent(inout) :: this
            real(kind=8),            intent(in)    :: a
        end subroutine smult_ifc
    end interface

    type, extends(stvec_abstract_t) :: stvec_t
        integer(kind=4)           :: N
        real(kind=8), allocatable :: p(:,:)
        contains
        !sum of two state vectors
        procedure             :: stvec_add_stvec
        !multiplication by scalar
        procedure             :: stvec_mult_scalar
        procedure, pass(this) :: scalar_mult_stvec
        !memory safe assign
        procedure, pass(lhs)  :: stvec_assign_stvec

        procedure             :: linear_comb
        procedure             :: smult
    end type stvec_t

contains
    subroutine stvec_assign_stvec_base(lhs, rhs)
        class(stvec_abstract_t), intent(inout) :: lhs
        class(stvec_abstract_t), intent(in)    :: rhs
        print *, "Error: assing not implemented for specific type"
        stop
    end subroutine stvec_assign_stvec_base

    function stvec_add_stvec(this, other) result(summed)
        class(stvec_t),          intent(in)  :: this
        class(stvec_abstract_t), intent(in)  :: other
        class(stvec_abstract_t), allocatable :: summed

        integer(kind=4) N

        select type(other)
        class is (stvec_t)
            N = this%N
            allocate(stvec_t :: summed)
            select type(summed)
            class is (stvec_t)
                allocate(summed%p(N,N))
                summed%p(1:N,1:N) = this%p(1:N,1:N)+other%p(1:N,1:N)
                summed%N = N
            end select
        class default
            print *, "error, stvec add type mismatch"
            stop
        end select
    end function stvec_add_stvec

    function scalar_mult_stvec(a, this) result(mult)
        class(stvec_t), intent(in)           :: this
        real(kind=8),   intent(in)           :: a
        class(stvec_abstract_t), allocatable :: mult

        integer(kind=4) N

        allocate(stvec_t :: mult)

        select type(mult)
        class is (stvec_t)
            N = this%N
            allocate(mult%p(N,N))
            mult%p(1:N,1:N) = a*this%p(1:N,1:N)
            mult%N = N
        end select
    end function scalar_mult_stvec

    function stvec_mult_scalar(this,a) result(mult)
        real(kind=8),   intent(in)           :: a
        class(stvec_t), intent(in)           :: this
        class(stvec_abstract_t), allocatable :: mult

        integer(kind=4) N

        allocate(stvec_t :: mult)

        select type(mult)
        class is (stvec_t)
            N = this%N
            allocate(mult%p(N,N))
            mult%p(1:N,1:N) = a*this%p(1:N,1:N)
            mult%N = N
        end select
    end function stvec_mult_scalar

    subroutine stvec_assign_stvec(lhs, rhs)
        class(stvec_t),           intent(inout) :: lhs
        class(stvec_abstract_t),  intent(in)    :: rhs

        integer(kind=4) N

        select type(rhs)
        class is(stvec_t)
            N = rhs%N
            if(.not.allocated(lhs%p)) allocate(lhs%p(N,N))
            lhs%p(1:N,1:N) = rhs%p(1:N,1:N)
            lhs%N = N
        class default
            !print *, "type mismatch in stvec assign"
            !stop
        end select

    end subroutine stvec_assign_stvec

    subroutine init_stvec(new_stvec, N, p)
        type(stvec_t),   intent(out) :: new_stvec
        integer(kind=4), intent(in)  :: N
        real(kind=8),    intent(in), &
                         optional    :: p(N, N)

        new_stvec%N = N
        allocate(new_stvec%p(N,N))
        if(present(p)) new_stvec%p(1:N,1:N) = p(1:N,1:N)
    end subroutine init_stvec

    subroutine linear_comb(this,other,a,b)
        class(stvec_t),          intent(inout) :: this
        class(stvec_abstract_t), intent(in)    :: other
        real(kind=8),            intent(in)    :: a, b

        integer(kind=4) N

        N = this%N

        select type (other)
        class is (stvec_t)
            this%p(1:N,1:N) = a*this%p(1:N,1:N)+b*other%p(1:N,1:N)
        class default
            print *, "type mismatch in stvec lincomb"
            stop
        end select
    end subroutine linear_comb
    subroutine smult(this,a)
        class(stvec_t), intent(inout) :: this
        real(kind=8),   intent(in)    :: a

        this%p = a*this%p
    end subroutine smult
end module stvec_mod
