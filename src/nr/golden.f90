    function golden(ax,bx,cx,func,tol,xmin)
    use nrtype
    implicit none
    real(dp), intent(in) :: ax,bx,cx,tol
    real(dp), intent(out) :: xmin
    real(dp) :: golden
    interface
        function func(x)
        use nrtype
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: func
        end function func
    end interface
    real(dp), parameter :: R=0.61803399_dp,C=1.0_dp-R
    real(dp) :: f1,f2,x0,x1,x2,x3
    x0=ax
    x3=cx
    if (abs(cx-bx) > abs(bx-ax)) then
        x1=bx
        x2=bx+C*(cx-bx)
    else
        x2=bx
        x1=bx-C*(bx-ax)
    end if
    f1=func(x1)
    f2=func(x2)
    do
        if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
        if (f2 < f1) then
            call shft3(x0,x1,x2,R*x2+C*x3)
            call shft2(f1,f2,func(x2))
        else
            call shft3(x3,x2,x1,R*x1+C*x0)
            call shft2(f2,f1,func(x1))
        end if
    end do
    if (f1 < f2) then
        golden=f1
        xmin=x1
    else
        golden=f2
        xmin=x2
    end if
    contains
!BL
    subroutine shft2(a,b,c)
    real(dp), intent(out) :: a
    real(dp), intent(inout) :: b
    real(dp), intent(in) :: c
    a=b
    b=c
    end subroutine shft2
!BL
    subroutine shft3(a,b,c,d)
    real(dp), intent(out) :: a
    real(dp), intent(inout) :: b,c
    real(dp), intent(in) :: d
    a=b
    b=c
    c=d
    end subroutine shft3
    end function golden
