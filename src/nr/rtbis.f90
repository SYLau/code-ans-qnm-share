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
    integer(I4B), parameter :: MAXIT=40
    integer(I4B) :: j
    real(dp) :: dx,f,fmid,xmid
    fmid=func(x2)
    f=func(x1)
    if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
    if (f < 0.0) then
        rtbis=x1
        dx=x2-x1
    else
        rtbis=x2
        dx=x1-x2
    end if
    do j=1,MAXIT
        dx=dx*0.5_dp
        xmid=rtbis+dx
        fmid=func(xmid)
        if (fmid <= 0.0) rtbis=xmid
        if (abs(dx) < xacc .or. fmid == 0.0) return
    end do
    call nrerror('rtbis: too many bisections')
end function rtbis
