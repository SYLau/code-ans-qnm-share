module nr_mod

    interface
        function golden(ax,bx,cx,func,tol,xmin)
        use nrtype
        real(dp), intent(in) :: ax,bx,cx,tol
        real(dp), intent(out) :: xmin
        real(dp) :: golden
        interface
            function func(x)
            use nrtype
            real(dp), intent(in) :: x
            real(dp) :: func
            end function func
        end interface
        end function golden
    end interface

    interface
        subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
            use nrtype
            real(dp), intent(inout) :: ax,bx
            real(dp), intent(out) :: cx,fa,fb,fc
            interface
                function func(x)
                    use nrtype
                    real(dp), intent(in) :: x
                    real(dp) :: func
                end function func
            end interface
        end subroutine mnbrak
    end interface

    interface
        subroutine polint(xa,ya,x,y,dy)
            use nrtype
            real(dp), dimension(:), intent(in) :: xa,ya
            real(dp), intent(in) :: x
            real(dp), intent(out) :: y,dy
        end subroutine polint
    end interface

    interface
        function rtsec(func,x1,x2,xacc)
            use nrtype
            real(dp), intent(in) :: x1,x2,xacc
            real(dp) :: rtsec
            interface
                function func(x)
                    use nrtype
                    real(dp), intent(in) :: x
                    real(dp) :: func
                end function func
            end interface
        end function rtsec
    end interface

    interface
        function rtbis(func,x1,x2,xacc)
        use nrtype; use nrutil, only : nrerror
        implicit none
        real(dp), intent(in) :: x1,x2,xacc
        real(dp) :: rtbis
        interface
            function func(x)
                use nrtype
                implicit none
                real(dp), intent(in) :: x
                real(dp) :: func
            end function func
        end interface
    end function rtbis
    end interface

    interface
        subroutine sort(arr)
        use nrtype
        real(dp), dimension(:), intent(inout) :: arr
        end subroutine sort
    end interface

    interface
        subroutine odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
            use nrtype; use nrutil, only : nrerror,reallocate
!            use ode_path
            implicit none
            real(dp), dimension(:), intent(inout) :: ystart
            real(dp), intent(in) :: x1,x2,eps,h1,hmin
            interface
                subroutine derivs(x,y,dydx)
                use nrtype
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
                end subroutine derivs
        !BL
                subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
                use nrtype
                implicit none
                real(dp), dimension(:), intent(inout) :: y
                real(dp), dimension(:), intent(in) :: dydx,yscal
                real(dp), intent(inout) :: x
                real(dp), intent(in) :: htry,eps
                real(dp), intent(out) :: hdid,hnext
                interface
                    subroutine derivs(x,y,dydx)
                    use nrtype
                    implicit none
                    real(dp), intent(in) :: x
                    real(dp), dimension(:), intent(in) :: y
                    real(dp), dimension(:), intent(out) :: dydx
                    end subroutine derivs
                end interface
                end subroutine rkqs
            end interface
        end subroutine odeint

        subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
            use nrtype; use nrutil, only : assert_eq,nrerror
            implicit none
            real(dp), dimension(:), intent(inout) :: y
            real(dp), dimension(:), intent(in) :: dydx,yscal
            real(dp), intent(inout) :: x
            real(dp), intent(in) :: htry,eps
            real(dp), intent(out) :: hdid,hnext
            interface
                subroutine derivs(x,y,dydx)
                use nrtype
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
                end subroutine derivs
            end interface
        end subroutine rkqs

        subroutine rkck(y,dydx,x,h,yout,yerr,derivs)
            use nrtype; use nrutil, only : assert_eq
            implicit none
            real(DP), dimension(:), intent(in) :: y,dydx
            real(DP), intent(in) :: x,h
            real(DP), dimension(:), intent(out) :: yout,yerr
            interface
                subroutine derivs(x,y,dydx)
                use nrtype
                implicit none
                real(DP), intent(in) :: x
                real(DP), dimension(:), intent(in) :: y
                real(DP), dimension(:), intent(out) :: dydx
                end subroutine derivs
            end interface
        end subroutine rkck
    end interface
endmodule nr_mod
