subroutine polint(xa,ya,x,y,dy)
    use nrtype; use nrutil, only : assert_eq,iminloc,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y,dy
    integer(I4B) :: m,n,ns
    real(dp), dimension(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
        den(1:n-m)=ho(1:n-m)-ho(1+m:n)
        if (any(den(1:n-m) == 0.0)) &
            call nrerror('polint: calculation failure')
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m)
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then
            dy=c(ns+1)
        else
            dy=d(ns)
            ns=ns-1
        end if
        y=y+dy
    end do
end subroutine polint
