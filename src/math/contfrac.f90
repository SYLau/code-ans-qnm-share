module contfrc_mod
    use type_mod,only:wp
    implicit none
    private
    public::mod_lentz

    interface
        function se_r(n) result(res)
            import wp
            integer,intent(in)::n
            real(wp)::res
        end function se_r
    end interface
    interface
        function se_c(n) result(res)
            import wp
            integer,intent(in)::n
            complex(wp)::res
        end function se_c
    end interface

    !=============================================================================================
    ! Using Modifed Lentz method to compute continued fraction in the form:
    ! f(j) = a1/(b1+) a2/(b2+) ... aj/bj
    !
    ! Ref: Numerical recipes
    !=============================================================================================
    type mod_lentz

!        real(wp)::eps = 1.e-15_wp
        real(wp)::eps = 1.e-6_wp
        real(wp)::tiny = 1.e-30_wp

        integer::itmax = 50

        procedure(se_r), pointer, nopass:: ar           => null()
        procedure(se_r), pointer, nopass:: br           => null()

        procedure(se_c), pointer, nopass:: ac           => null()
        procedure(se_c), pointer, nopass:: bc           => null()

        contains

        procedure, pass :: r            => mod_lentz_r
        procedure, pass :: c            => mod_lentz_c

    end type mod_lentz

    contains

    function mod_lentz_r(this) result(f)
        class(mod_lentz)::this
        real(wp)::f
        real(wp)::C,D
        integer::j

        if (.not.associated(this%ar) .or. .not.associated(this%br)) error stop 'err: mod_lentz_r series not associated'

        f = this%tiny
        C = f
        D = 0
        do j=1,this%itmax
            D = this%br(j) + this%ar(j)*D
            if (abs(D) <= epsilon(abs(D))) D = this%tiny
            D = 1/D

            C = this%br(j) + this%ar(j)/C
            if (abs(C) <= epsilon(abs(C))) C = this%tiny

            f = f*C*D

            if (abs(C*D - 1) < this%eps) return
        end do

!        error stop 'err: mod_lentz_r exceeds max iteration'

    end function mod_lentz_r

    function mod_lentz_c(this) result(f)
        class(mod_lentz)::this
        complex(wp)::f
        complex(wp)::C,D
        integer::j

        if (.not.associated(this%ac) .or. .not.associated(this%bc)) error stop 'err: mod_lentz_c series not associated'

        f = this%tiny
        C = f
        D = 0
        do j=1,this%itmax
            D = this%bc(j) + this%ac(j)*D
            if (abs(D) <= epsilon(abs(D))) D = this%tiny
            D = 1/D

            C = this%bc(j) + this%ac(j)/C
            if (abs(C) <= epsilon(abs(C))) C = this%tiny

            f = f*C*D

            if (abs(C*D - 1) < this%eps) return
        end do

!        error stop 'err: mod_lentz_c exceeds max iteration'

    end function mod_lentz_c

end module contfrc_mod
