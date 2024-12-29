module leaver_mod
    use type_mod,only:wp
    implicit none
    private
    public::RW_Zerilli
    public::leaver_solve


    contains

    function RW_Zerilli(l,r,m,wc,Z) result(Q)
        real(wp),intent(in)::l,r,m
        complex(wp),intent(in)::wc
        complex(wp),intent(in)::Z(2)
        complex(wp)::Q(2)
        complex(wp)::im
        real(wp)::kk,bb,ff,lam

        im=cmplx(0._wp,1._wp,kind=wp)

        lam=l*(l+1)/2-1
        kk=4*lam*(lam+1)
        bb=6*m
        ff=(r-2*m)/(2*r**2*(lam*r+3*m))

        Q(1) = ((kk + 2*bb**2*ff)*Z(1) - 2*bb*(1-2*m/r)*Z(2))/(kk-2*im*wc*bb)
        Q(2) = ((kk+2*im*wc*bb)*Z(1) - (kk+2*bb**2*ff)*Q(1))/(2*bb*(1-2*m/r))
    end function RW_Zerilli

    function leaver_solve(l,r,m,wc,Q) result(soln)
        use contfrc_mod,only: mod_lentz
        real(wp),intent(in)::l,r,m
        complex(wp),intent(in)::wc
        complex(wp),intent(in)::Q(2)
        complex(wp)::soln
        complex(wp)::im
        complex(wp)::beta0
        type(mod_lentz)::cf

        im=cmplx(0._wp,1._wp,kind=wp)

        beta0=r/Q(1)*(Q(2)+im*wc*r/(r-2*m)*Q(1))

        cf%itmax = 15

        cf%ac => ac
        cf%bc => bc
        soln = beta0 + cf%c()

        contains

        function ac(n) result(res)
            integer,intent(in) :: n
            complex(wp)::res

            res = -alpha(n-1)*ga(n)
        end function ac

        function bc(n) result(res)
            integer,intent(in) :: n
            complex(wp)::res

            res = beta(n)
        end function bc

        function alpha(n) result(res)
            integer,intent(in) :: n
            complex(wp)::res

            if (n<0) error stop 'err: leaver_solve alpha n < 0'
            if (n == 0) then
                res = -1._wp
            else
                res = (1 - 2*m/r)*n*(n+1)
            endif
        end function alpha

        recursive function beta(n) result(res)
            integer,intent(in) :: n
            complex(wp)::res
            real(wp)::delta

            if (n<=0) error stop 'err: leaver_solve beta n <= 0'
            if (n == 1) then
                res = -2*(im*wc*r+(1-3*m/r))
            else
                delta = 2*m/r*(n-3)*(n+1)
                res = -2*(im*wc*r+(1-3*m/r)*n)*n-alpha(n-1)*delta/ga(n-1)
            endif
        end function beta

        recursive function ga(n) result(res)
            integer,intent(in) :: n
            complex(wp)::res
            real(wp)::delta

            if (n<=0) error stop 'err: leaver_solve ga n <= 0'
            if (n == 1) then
                res = 6*m/r-l*(l+1)
            else
                delta = 2*m/r*(n-3)*(n+1)
                res = (1-6*m/r)*n*(n-1) + 6*m/r -l*(l+1) - beta(n-1)*delta/ga(n-1)
            endif
        end function ga


        end function leaver_solve
end module leaver_mod
