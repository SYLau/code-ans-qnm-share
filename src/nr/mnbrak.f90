subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
    use nrtype; use nrutil, only : swap
    implicit none
    real(dp), intent(inout) :: ax,bx
    real(dp), intent(out) :: cx,fa,fb,fc
    interface
        function func(x)
            use nrtype
            implicit none
            real(dp), intent(in) :: x
            real(dp) :: func
        end function func
    end interface
    real(dp), parameter :: GOLD=1.618034_dp,GLIMIT=100.0_dp,tiny=1.0e-20_dp
    real(dp) :: fu,q,r,u,ulim
    fa=func(ax)
    fb=func(bx)
    if (fb > fa) then
        call swap(ax,bx)
        call swap(fa,fb)
    end if
    cx=bx+GOLD*(bx-ax)
    fc=func(cx)
    do
        if (fb < fc) return
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),tiny),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if ((bx-u)*(u-cx) > 0.0) then
            fu=func(u)
            if (fu < fc) then
                ax=bx
                fa=fb
                bx=u
                fb=fu
                return
            else if (fu > fb) then
                cx=u
                fc=fu
                return
            end if
            u=cx+GOLD*(cx-bx)
            fu=func(u)
        else if ((cx-u)*(u-ulim) > 0.0) then
            fu=func(u)
            if (fu < fc) then
                bx=cx
                cx=u
                u=cx+GOLD*(cx-bx)
                call shft(fb,fc,fu,func(u))
            end if
        else if ((u-ulim)*(ulim-cx) >= 0.0) then
            u=ulim
            fu=func(u)
        else
            u=cx+GOLD*(cx-bx)
            fu=func(u)
        end if
        call shft(ax,bx,cx,u)
        call shft(fa,fb,fc,fu)
    end do
contains
    !BL
    subroutine shft(a,b,c,d)
        real(dp), intent(out) :: a
        real(dp), intent(inout) :: b,c
        real(dp), intent(in) :: d
        a=b
        b=c
        c=d
    end subroutine shft
end subroutine mnbrak
