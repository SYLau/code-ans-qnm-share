module zerilli_mod
    implicit none
    private
    public:: zerilli_eq,zerilli_initial,zerilli_asym_coef

contains
    subroutine zerilli_eq(l,wc2,m0,t,x,f)
        real(8), intent(in):: l,m0,t
        complex(8),intent(in)::wc2,x(:)
        complex(8), intent(out):: f(:)
        real(8)::nn,b1,b2,b3,b4,b5
        complex(8)::b(2,2)
        real(8)::r,drt,Vz

        nn=(l-1.d0)*(l+2.d0)/2.d0
        r=t
        drt=1.d0/(1.d0-2.d0*m0/r)
        b1=(1.d0-2.d0*m0/r)/(r**3 * (nn*r + 3.d0 * m0)**2)
        b2=2.d0*nn**2*(nn+1.d0)*r**3
        b3=6.d0*nn**2*m0*r**2
        b4=18.d0*nn*m0**2*r
        b5=18.d0*m0**3
        Vz=b1*(b2+b3+b4+b5)
        b(1,1)=0.d0
        b(1,2)=(1.d0,0.d0)
        b(2,1) = (Vz-wc2)*drt**2
        b(2,2) = -2.d0*m0/r**2*drt

        f(1)=dot_product(conjg(b(1,:)),x(:))
        f(2)=dot_product(conjg(b(2,:)),x(:))
    end subroutine zerilli_eq

    function zerilli_initial(m,r,H0,K,l,w)  !modified from pes03_bc_Oui; Initial condition at r=R0 for Zerilli eq from LD2
        use type_mod,only:c
        real(8),intent(in)::m,r,l
        complex(8),intent(in)::H0,K,w
        complex(8)::zerilli_initial(1:2)
        real(8) :: n, gg, ll, hh, kk, det, tcf
        complex(8) :: a,b, wc2

        tcf = (1.d0 - 2.d0*m/r)**(-1.d0)
        wc2=(w/c)**2

        n = (l-1.d0)*(l+2.d0)/2.d0
        a = -(n*r + 3.d0*m)/(wc2*r**2.d0 - (n+1.d0)*m/r)
        b = (n*r*(r-2.d0*m) - wc2*r**4.d0 + m*(r-3.d0*m))/(r-2.d0*m)/(wc2*r**2.d0 - (n+1.d0)*m/r) ! this mistake caused ~15% error in f-mode
        gg = (n*(n+1.d0)*r**2 + 3.d0*n*m*r + 6.d0*m**2)/r**2.d0/(n*r + 3.d0 * m)
        ll = 1.d0
        hh = (-n*r**2.d0 + 3.d0*n*m*r + 3.d0*m**2.d0)/(r-2.d0*m)/(n*r + 3.d0 * m)
        kk = -r**2.d0/(r - 2.d0*m)
        det = gg * kk - hh * ll
        zerilli_initial(1) = (-a*ll * H0 + (kk - b * ll) * K )/det
        zerilli_initial(2) = (gg *a * H0 + (-hh + b * gg) * K )/det*tcf
    end function

    subroutine zerilli_asym_coef(l,m0,r,w,x,be,ga)
        implicit none
        real(8),intent(in)::l,m0,r
        complex(8),intent(in)::w,x(1:2)
        complex(8),intent(out)::be,ga
        real(8) :: nn, rt, drt
        complex(8) ::  im, det, alp(0:2), calp(0:2)
        complex(8) :: z_asym(1:2, 1:2)

        rt = r + 2.d0*m0* log(r/(2.d0*m0) - 1.d0)
        im = (0.d0, 1.d0)
        nn = (l-1.d0)*(l+2.d0)/2.d0
        drt = 1.d0/(1.d0-2.d0*m0/r)

        alp(0) = 1.d0
        alp(1) = -im*(nn + 1.d0)/w *alp(0)
!        alp(2) = -0.5d0/w**2*(nn*(nn+1.d0) - 1.5d0*im*m0*w*(1.d0 + 2.d0/nn))*alp(0)
        !This expression contains a typo and is not reported until Lu and Suen 2011 https://iopscience.iop.org/article/10.1088/1674-1056/20/4/040401
        !Changed to:
        alp(2) = -1.d0/w**2*(0.5d0*nn*(nn+1.d0) - 1.5d0*im*m0*w*(1.d0 + 2.d0/nn))*alp(0)
        calp = conjg(alp)

        z_asym(1,1) = exp(-im*w*rt)*(alp(0) + alp(1)*r**(-1.d0) + alp(2)*r**(-2.d0))
        z_asym(1,2) = exp(im*w*rt)*(calp(0) + calp(1)*r**(-1.d0) + calp(2)*r**(-2.d0))
        z_asym(2,1) = -im*w*z_asym(1,1) + exp(-im*w*rt)*(-alp(1)*r**(-2.d0)-2.d0*alp(2)*r**(-3.d0))/drt
        z_asym(2,2) = im*w*z_asym(1,2) + exp(im*w*rt)*(-calp(1)*r**(-2.d0)-2.d0*calp(2)*r**(-3.d0))/drt

        det = z_asym(1,1)*z_asym(2,2) - z_asym(1,2)*z_asym(2,1)

        be = (x(1)*z_asym(2,2) - z_asym(1,2)*x(2)/drt )/det
        ga = (z_asym(1,1)*x(2)/drt - x(1)*z_asym(2,1) )/det
    endsubroutine zerilli_asym_coef

end module
