function rtsec(func,x1,x2,xacc)
    use nrtype; use nrutil, only : nrerror,swap
    implicit none
    real(dp), intent(in) :: x1,x2,xacc
    real(dp) :: rtsec
    interface
        function func(x)
            use nrtype
            implicit none
            real(dp), intent(in) :: x
            real(dp) :: func
        end function func
    end interface
    integer(I4B), parameter :: MAXIT=30
    integer(I4B) :: j
    real(dp) :: dx,f,fl,xl
    fl=func(x1)
    f=func(x2)
    if (abs(fl) < abs(f)) then
        rtsec=x1
        xl=x2
        call swap(fl,f)
    else
        xl=x1
        rtsec=x2
    end if
    do j=1,MAXIT
        dx=(xl-rtsec)*f/(f-fl)
        xl=rtsec
        fl=f
        rtsec=rtsec+dx
        f=func(rtsec)
        if (abs(dx) < xacc .or. f == 0.0) return
    end do
    call nrerror('rtsec: exceed maximum iterations')
end function rtsec
