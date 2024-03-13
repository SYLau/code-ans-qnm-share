module muller_mod
    implicit none
    private
    public::rtmuller

    !precision
    integer, parameter::dp=kind(1.d0)
    integer, parameter::qp=kind(1.q0)
    integer, parameter::wp=dp

contains
    function rtmuller(func, x1_in, x2_in, x3_in, eps, quit_if_error, ostat, maxit)
        implicit none
        complex(wp),intent(in)::x1_in, x2_in, x3_in
        real(wp),intent(in)::eps
        logical, intent(in), optional::quit_if_error
        integer,intent(out),optional::ostat
        integer,intent(in),optional::maxit
        complex(wp) :: x1, x2, x3, f1, f2, f3, xtemp!, dx1, dx2
        complex(wp) :: c, q, a, b, discrim, denom
        real(wp) :: x2_bar
        complex(wp) :: rtmuller
        integer :: n
        integer :: JMAX=100
        interface
            function func(x)
                import wp
                implicit none
                complex(wp), intent(in) :: x
                complex(wp) :: func
            end function func
        end interface

        if (present(maxit)) JMAX=maxit

        n = 0
        x1 = x1_in
        x2 = x2_in
        x3 = x3_in

        !dx1 = x2 - x1
        !dx2 = x3 - x2
        x2_bar = (abs(x3) + abs(x2))/2._wp
        if (x1 == x2) then
            write(*,*) "err: rtmuller dx=0"
            stop
        end if

        do
            n = n + 1

            do while ((x3-x2)*(x2-x1)*(x1-x3) == 0._wp)
                if (x3-x2 == 0._wp) x3 = x3*1.00001_wp
                if (x2-x1 == 0._wp) x2 = x2*1.00002_wp
                if (x1-x3 == 0._wp) x1 = x1*1.00003_wp
            enddo

            f1 = func(x1)
            f2 = func(x2)
            f3 = func(x3)

            if (f1 == 0._wp) then
                xtemp = x1
                exit
            elseif (f2 == 0._wp) then
                xtemp = x2
                exit
            elseif (f3 == 0._wp) then
                xtemp = x3
                exit
            endif

            q=(x3-x2)/(x2-x1)
            a = q*f3-q*(1._wp-q)*f2+q**2*f1
            b = (2._wp*q+1._wp)*f3-(1._wp+q)**2*f2+q**2*f1
            c = (1._wp+q)*f3

            discrim = sqrt(b*b - 4._wp * a * c)

            IF ( abs(b+discrim) > abs(b-discrim)) THEN
                denom = b + discrim
            ELSE
                denom = b - discrim
            END IF

            xtemp = x3 - (x3-x2)*(2._wp*c) / denom

            !dx1 = x2 - x1
            !dx2 = x3 - x2
            x2_bar = (abs(x3) + abs(x2))/2._wp
            x1 = x2
            x2 = x3
            x3 = xtemp
!print*,abs(x3-x2),eps* abs(x2_bar),eps* abs(x2_bar)/abs(x3-x2)
            if (abs(x3-x2) < eps* abs(x2_bar)) exit
            if (n >= JMAX) then
                write(*,*) "err: rtmuller too many steps"
                if (present(quit_if_error)) then
                    if (quit_if_error) then
                        stop
                    else
                        if (present(ostat)) ostat=1
                        return
                    end if
                end if
!                pause "err: rtmuller too many steps"
!                return
            endif
        enddo

        rtmuller = xtemp
        if (present(ostat)) ostat=0

    endfunction rtmuller

!    function rtmuller_backup(func, x1_in, x2_in, x3_in, eps)
!        implicit none
!        complex(wp),intent(in)::x1_in, x2_in, x3_in
!        real(wp),intent(in)::eps
!        complex(wp) :: x1, x2, x3, f1, f2, f3, xtemp, dx1, dx2
!        complex(wp) :: c, t, u, v, a, b, discrim, denom
!        real(wp) :: x2_bar
!        complex(wp) :: rtmuller_backup
!        integer :: n
!        integer, parameter :: JMAX=30
!        interface
!            function func(x)
!                implicit none
!                complex(wp), intent(in) :: x
!                complex(wp) :: func
!            end function func
!        end interface
!
!        n = 0
!        x1 = x1_in
!        x2 = x2_in
!        x3 = x3_in
!
!        dx1 = x2 - x1
!        dx2 = x3 - x2
!        x2_bar = (abs(x3) + abs(x2))/2._wp
!        if (dx1 == (0._wp,0._wp) .or. dx2 == (0._wp,0._wp)) pause "err: rtmuller dx=0"
!
!        do
!            n = n + 1
!
!            do while ((x3-x2)*(x2-x1)*(x1-x3) == 0._wp)
!                if (x3-x2 == 0._wp) x3 = x3*1.00001_wp
!                if (x2-x1 == 0._wp) x2 = x2*1.00002_wp
!                if (x1-x3 == 0._wp) x1 = x1*1.00003_wp
!            enddo
!
!            f1 = func(x1)
!            f2 = func(x2)
!            f3 = func(x3)
!
!            if (f1 == 0._wp) then
!                xtemp = x1
!                exit
!            elseif (f2 == 0._wp) then
!                xtemp = x2
!                exit
!            elseif (f3 == 0._wp) then
!                xtemp = x3
!                exit
!            endif
!
!            c = f3
!            t = (f3-f2)/(x3-x2)
!            u = (f2-f1)/(x2-x1)
!            v = (f1-f3)/(x1-x3)
!            a = (v-u)/(x3-x2)
!            b = t + v - u
!
!            discrim = sqrt(b*b - 4._wp * a * c)
!
!            IF ( abs(b+discrim) > abs(b-discrim)) THEN
!                denom = b + discrim
!            ELSE
!                denom = b - discrim
!            END IF
!
!            xtemp = x3 - (2._wp*c) / denom
!
!            dx1 = x2 - x1
!            dx2 = x3 - x2
!            x2_bar = (abs(x3) + abs(x2))/2._wp
!            x1 = x2
!            x2 = x3
!            x3 = xtemp
!            if (abs(dx2) < eps* abs(x2_bar)) exit
!            if (n >= JMAX) then
!                pause "err: rtmuller too many steps"
!                return
!            endif
!        enddo
!        rtmuller_backup = xtemp
!    endfunction rtmuller_backup

end module muller_mod

