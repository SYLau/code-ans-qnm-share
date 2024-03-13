module ld_eq_aniso_mod
    use type_mod,only:wp
    !The module contains the equations from LD formalism
    !all variables in geometrized unit
    implicit none
    private
    public::def_sigma
    public::perturb_par
    public::ld_eq,get_H0_V_Z,ld_r0_reg

    interface
        function sigma_func(p,rho,m,r) result(res)
            import wp
            real(wp),intent(in)::p,rho,m,r
            real(wp)::res
        end function sigma_func

        function dsigma_func(p,rho,m,ga,A,r) result(res)
            import wp
            real(wp),intent(in)::p,rho,m,ga,A,r
            real(wp)::res
        end function dsigma_func

        function sigma_coef(p,rho,ga,A) result(res)
            import wp
            real(wp),intent(in)::p,rho,ga,A
            real(wp)::res
        end function sigma_coef
    end interface

    type perturb_par
        real(wp)::l,p,rho,m,nu,ga,eLam,dLam,dnu,eNu
        real(wp)::A
        complex(wp)::wc2
    end type perturb_par

    procedure (sigma_func), pointer:: get_sigma             => null()
    procedure (dsigma_func), pointer:: get_dsigma_p          => null()
    procedure (dsigma_func), pointer:: get_dsigma_rho        => null()
    procedure (dsigma_func), pointer:: get_dsigma_mu         => null()

    procedure (dsigma_func), pointer:: get_ds                => null()

    procedure (sigma_coef), pointer:: get_sigma_2           => null()

    logical:: has_sigma = .false.


!    real(wp)::temp_par=0

contains

    subroutine def_sigma(s,ds_p,ds_rho,ds_mu,ds,s2)
        interface
            function s(p,rho,m,r) result(res)
                import wp
                real(wp),intent(in)::p,rho,m,r
                real(wp)::res
            end function s
            function ds_p(p,rho,m,ga,A,r) result(res)
                import wp
                real(wp),intent(in)::p,rho,m,ga,A,r
                real(wp)::res
            end function ds_p
            function ds_rho(p,rho,m,ga,A,r) result(res)
                import wp
                real(wp),intent(in)::p,rho,m,ga,A,r
                real(wp)::res
            end function ds_rho
            function ds_mu(p,rho,m,ga,A,r) result(res)
                import wp
                real(wp),intent(in)::p,rho,m,ga,A,r
                real(wp)::res
            end function ds_mu
            function ds(p,rho,m,ga,A,r) result(res)
                import wp
                real(wp),intent(in)::p,rho,m,ga,A,r
                real(wp)::res
            end function ds
            function s2(p,rho,ga,A) result(res)
                import wp
                real(wp),intent(in)::p,rho,ga,A
                real(wp)::res
            end function s2
        end interface

        get_sigma => s
        get_dsigma_p => ds_p
        get_dsigma_rho => ds_rho
        get_dsigma_mu => ds_mu

        get_ds => ds

        get_sigma_2 => s2

        has_sigma = .true.
    end subroutine def_sigma

    !=======================================================================
    ! LD equations modified for anisotropic case
    !=======================================================================
    subroutine ld_eq(par,t,xc,fc)
        use type_mod,only:pi
        type(perturb_par),intent(in)::par
        real(wp),intent(in)::t
        complex(wp),intent(in)::xc(4)
        complex(wp), intent(out)::fc(4)
        real(wp)::l,p,rho,m,nu,ga,eLam,dLam,dnu,eNu
        real(wp)::A
        complex(wp)::wc2
        real(wp)::s
        real(wp)::ds_dr
        real(wp)::sb
        real(wp)::r
        complex(wp)::H1,K,W,X
        complex(wp)::H0,V,Z

        l=par%l; p=par%p; rho=par%rho; m=par%m; nu=par%nu; ga=par%ga
        A=par%A
        eLam=par%eLam; dLam=par%dLam; dNu=par%dnu; eNu=par%eNu
        wc2=par%wc2

        r=t

        s=get_sigma(p,rho,m,r)
        ds_dr= get_ds(p,rho,m,ga,A,r)

        sb=s/(rho+p)

        H1=xc(1); K=xc(2); W=xc(3); X=xc(4)
        call get_H0_V_Z(par,t,xc,H0,V,Z)

        H1_eq: block
            complex(wp)::coef_H1,coef_K,coef_H0,coef_V

            coef_H1= pi*4*(rho-p)*eLam*r - 2*m/r**2*eLam - (l+1)/r
            coef_K= eLam/r
            coef_H0= eLam/r
            coef_V= -16*pi*(rho+p)/r*eLam*(1-sb)

            fc(1) = coef_H1*H1 + coef_K*K +coef_H0*H0 + coef_V*V
        end block H1_eq

        K_eq: block
            complex(wp)::coef_H1,coef_K,coef_W,coef_H0

            coef_H1=l*(l+1)/2/r
            coef_K=dNu/2-(l+1)/r
            coef_W=-pi*8*(rho+p)*sqrt(eLam)/r
            coef_H0=1/r

            fc(2) = coef_H1*H1 + coef_K*K +coef_W*W +coef_H0*H0
        end block K_eq

        W_eq: block
            complex(wp)::coef_K,coef_W,coef_X,coef_H0,coef_V

            coef_K= r*sqrt(eLam)*(1-sb)
            coef_W= -(l+1)/r +2*sb/r
            coef_X= r*sqrt(eLam/eNu)/ga/p
            coef_H0= r*sqrt(eLam)/2
            coef_V= -l*(l+1)/r*sqrt(eLam)*(1-sb)

            fc(3) = coef_K*K +coef_W*W  +coef_X*X +coef_H0*H0 +coef_V*V
        end block W_eq

        X_eq: block
            complex(wp)::coef_H1,coef_K,coef_W,coef_X,coef_H0,coef_V,coef_Z
            real(wp)::dp
            real(wp)::t1,t2

            dp= -(rho+p)/2*dNu -2*s/r
            t1= (rho+p)/2*sqrt(eNu)
            t2= eLam/r**4*(7*m**2 -4*r*m -8*pi*r**3*m*rho -16*pi**2*r**6*p**2 -pi*4*(p-rho)*r**4 - 8*pi*s/eLam*r**4) &
             - ( (6/r**2-2*dNu/r)*sb -2/r**2*r*ds_dr/(rho+p) - 4/r**2*sb**2)/eLam
            coef_H1= t1*(r*wc2/eNu + l*(l+1)/2/r*(1-2*sb))
            coef_K= t1*( (1.5_wp-2*sb)*dNu - (1-6*sb)/r -4*sb**2/r)
            coef_W= -t1*2/r*sqrt(eLam)*(4*pi*(rho+p) +wc2/eNu - t2)
            coef_X=-1/r*(l-2*(rho+p)/ga/p*sb)
            coef_H0=t1*(1/r-dNu/2)
            coef_V=l*(l+1)*sqrt(eNu)*dp/r**2*(1-sb)
            coef_Z=2*sqrt(eNu)/r

            fc(4) = coef_H1*H1 + coef_K*K +coef_W*W  +coef_X*X +coef_H0*H0 +coef_V*V +coef_Z*Z
        end block X_eq
    end subroutine ld_eq

    subroutine get_H0_V_Z(par,t,xc,H0,V,Z) !Solving LD2 Eq (5-6) for H0 and V in terms of {H1,K,W,X}
        use type_mod,only:pi
        type(perturb_par),intent(in)::par
        real(wp), intent(in)::t
        complex(wp),intent(in)::xc(4)
        complex(wp),intent(out)::H0,V,Z
        real(wp)::l,p,rho,m,nu,ga,eLam,dLam,dnu,eNu
        real(wp)::A
        complex(wp)::wc2
        real(wp)::r
        real(wp)::s
        real(wp)::sb
        complex(wp)::H1,K,W,X

        l=par%l; p=par%p; rho=par%rho; m=par%m; nu=par%nu; ga=par%ga
        A=par%A
        eLam=par%eLam; dLam=par%dLam; dNu=par%dnu; eNu=par%eNu
        wc2=par%wc2

        r=t

        s=get_sigma(p,rho,m,r)

        sb= s/(rho+p)

        H1=xc(1); K=xc(2); W=xc(3); X=xc(4)

        H0_eq: block
            complex(wp):: coef_H1, coef_K, coef_W, coef_X
            real(wp)::t1

            t1=3*m+(l+2)*(l-1)/2*r+4*pi*r**3*p
            coef_H1= -(l*(l+1)/2*(m+4*pi*r**3*p) -wc2/eNu *r**3/eLam)/t1
            coef_K= ((l+2)*(l-1)/2*r -wc2/eNu*r**3 -eLam/r*(m+4*pi*r**3*p)*(3*m-r+4*pi*r**3*p) )/t1
            coef_W= -16*pi*r/sqrt(eLam)*(rho+p)*sb/t1
            coef_X= 8*pi*r**3/sqrt(eNu)/t1

            H0 = coef_H1*H1 + coef_K*K + coef_W*W + coef_X*X
        end block H0_eq

        Z_eq: block
            complex(wp):: coef_W, coef_X, coef_H0
            real(wp)::dp
            real(wp)::cs2
            real(wp)::dsdp,dsdrho,dsdmu

            cs2 = p/(rho+p)*ga
            dp = -(rho+p)/2*dNu - 2*s/r
            dsdp = get_dsigma_p(p,rho,m,ga,A,r)
            dsdrho = get_dsigma_rho(p,rho,m,ga,A,r)
            dsdmu = get_dsigma_mu(p,rho,m,ga,A,r)

            coef_W= -dsdp*dp/r/sqrt(eLam) -dsdrho*((rho+p)*A/r+dp/r/cs2)/sqrt(eLam)
            coef_X= -dsdp/sqrt(eNu) -dsdrho/cs2/sqrt(eNu)
            coef_H0= -dsdmu/eLam

!            coef_W= -dsdp*dp/r/sqrt(eLam)
!            coef_X= -dsdp/sqrt(eNu)
!            coef_H0= -dsdmu/eLam

            Z = coef_W*W + coef_X*X + coef_H0*H0
        end block Z_eq

        V_eq: block
            complex(wp):: coef_W, coef_X, coef_H0, coef_Z
            real(wp)::dp
            complex(wp)::t1

            dp = -(rho+p)/2*dNu - 2*s/r
            t1 = wc2/sqrt(eNu)*(rho+p)*(1-sb)

            coef_W=dp/r*sqrt(eNu/eLam)/t1
            coef_X=1/t1
            coef_H0=-(rho+p)/2*sqrt(eNu)/t1
            coef_Z=sqrt(eNu)/t1

            V = coef_W*W + coef_X*X + coef_H0*H0 + coef_Z*Z
        end block V_eq

    end subroutine

    function ld_r0_reg(isol,P0,rho0,Nu0,gamma0,r,l,wc) !modified from pes03_bc_Co_Coef; Regular sol at r=0 from LD2; isol=1,2 represents 2 linearly indep. sols.
        use type_mod,only:pi
!        use linsys_mod,only:LinSys_Cramer
        integer,intent(in)::isol
        real(wp),intent(in)::P0, rho0, Nu0, gamma0, r,l
        complex(wp),intent(in)::wc
        complex(wp)::ld_r0_reg(4)
        real(wp) :: l1, l2, pi4_3, pi4_5, eNu0
        real(wp) :: P2, rho2, Nu2, P4, Nu4
        real(wp) :: s2
!real(wp) :: ds_dmu
!        complex(wp) :: Q0, Q1
        complex(wp) :: H10, K0, W0, X0
!complex(wp) :: Z0
        complex(wp) :: b1,b2,b3,b4, a(2,4), wc2b
        complex(wp) :: cs2,Z2
        complex(wp) :: H12, K2, W2, X2

        s2 = get_sigma_2(P0, rho0, gamma0, 0.d0)
!ds_dmu=get_dsigma_mu(p0,rho0,pi*4/3*rho0*r**3,r)
        pi4_3 = 4.d0/3*pi
        pi4_5 = 4.d0/5*pi
        l1 = l*(l+1)
        l2 = (l+2)*(l-1)

        eNu0 = exp(Nu0)**2
        wc2b = wc**2*eNu0**(-1.)

        P2 = -pi4_3*(P0 + rho0)*(3*P0 + rho0)
        rho2 = P2*(P0 + rho0)/gamma0/P0
        Nu2 = pi4_3*2*(3*P0 + rho0)

        P4 = -pi4_5/2*(P0+rho0)*(5.d0*P2+rho2)-pi4_3/2*(P2+rho2)*(3*P0+rho0) &
            -32._wp/9*pi**2*rho0*(P0+rho0)*(3*P0+rho0)
        Nu4 = pi4_5 * (5*P2 + rho2) + 4*pi4_3**2 * rho0 * (3*P0 + rho0)

        if (isol == 1) K0 = cmplx(P0+rho0, 0.,kind=8)
        if (isol == 2) K0 = -cmplx(P0+rho0, 0.,kind=8)

!Z0 = ds_dmu*K0
!Z0 = 0

        W0 = (1._wp, 0._wp)
        H10 = 1._wp/l1*(2*l*K0 + 16*pi*(P0+rho0)*W0)
        b1 = (P0 + rho0)*eNu0**(0.5_wp)
        b2 = 0.5_wp
        b3 = pi*4/3*(3*P0+rho0)-wc2b/l + 2*s2/(rho0+P0)
!        X0 = b1 *(b2 * K0 + b3 * W0) - eNu0**(0.5_wp)*Z0
        X0 = b1 *(b2 * K0 + b3 * W0)

        a(1,1) = H10
        a(1,2) = K0
        a(1,3) = W0
        a(1,4) = X0

        cs2 = gamma0*P0/(P0 + rho0)

Z2 = 0

        H12 = (-wc2b/14+pi/63*(58*rho0-114*p0)+pi*8/7*s2/wc2b*(1-1._wp/cs2))*K0 -pi*48/7/wc2b*Z2 &
        +(-pi*2/7*wc2b+16*pi**2/63*(19*rho0-15*p0)-16*pi**2/21/cs2*(rho0+3*p0))*(rho0+p0)*W0

!        Q0 = 4.d0/l2 *(pi4*2.d0 * eNu0**(-0.5d0) * X0 - (pi4_3*2.d0 * rho0 + wc2b)*K0 &
!            -(pi4_3*0.5d0*l1*(3.d0*P0 + rho0) - wc2b)*H10)
!        Q1 = 2.d0/l1 * (eNu0**(-0.5d0)/gamma0/P0* X0 + 1.5d0* K0 + pi4_3*(l + 1.d0)*rho0 *W0)
!
!        U(1,1) = 0.d0
!        U(1,2) = -(P0 + rho0)/4.d0
!        U(1,3) = 0.5d0 * (P2 + (P0 + rho0)*wc2b*(l+3.d0)/l1 )
!        U(1,4) = 0.5d0* eNu0**(-0.5d0)
!
!        U(2,1) = - l1/4.d0
!        U(2,2) = (l + 2.d0)/2.d0
!        U(2,3) = pi4 * (P0 + rho0)
!        U(2,4) = 0.d0
!
!        U(3,1) = (l + 3.d0)/2.d0
!        U(3,2) = -1.d0
!        U(3,3) = -pi4*2.d0 * (P0 + rho0)*(l+3.d0)/l1
!        U(3,4) = 0.d0
!
!        U(4,1) = -l1 * (P0 + rho0) * eNu0**(0.5d0) /8.d0
!        U(4,2) = 0.d0
!        U(4,3) = -(P0 + rho0) * eNu0**(0.5d0) * ((l + 2.d0)/4.d0 * Nu2 - pi4/2.d0 * (P0 + rho0) - wc2b/2.d0)
!        U(4,4) = (l + 2.d0)/2.d0
!
!        b1 = Nu2*eNu0**(-0.5d0)/4.d0 * X0
!        b2 = (rho2 + P2)/4.d0 * K0
!        b3 = (P0 + rho0)/4.d0* Q0 + wc2b*(P0 + rho0)/2.d0 * Q1
!        b4 = -(P4 - pi4_3*rho0*P2 + wc2b/2.d0/l*(rho2 + P2 - (P0 + rho0)*Nu2 ) )* W0
!        Z(1) = b1+b2+b3+b4
!        b1 = pi4_3 * (3.d0*P0 + rho0) * K0
!        b2 = 0.5d0 * Q0
!        b3 = -pi4 * (rho2 + P2 + pi4_3*2.d0 * rho0 * (P0 * rho0))*W0
!        Z(2) = b1+b2+b3
!        b1 = pi4*((2.d0*l + 3.d0)*rho0/3.d0 - P0)*H10
!        b2 = pi4*2.d0/l * (P2 + rho2)*W0
!        b3 = -pi4*2.d0 * (P0 + rho0)* Q1 + 0.5d0 * Q0
!        Z(3) = b1+b2+b3
!        b1 = 0.5d0*(rho2 + P2 + 0.5d0*(P0 + rho0)*Nu2)*l/(P0 + rho0) * X0
!        b2 = (P0 + rho0)*eNu0**(0.5d0) * (0.5d0*Nu2* K0 + 0.25d0* Q0 + 0.5d0*wc2b* H10 -0.25d0*l1*Nu2* Q1)
!        b3 = (P0 + rho0)*eNu0**(0.5d0) * (0.5d0*(l + 1.d0)*Nu4 - pi4/2.d0*(P2 + rho2) - pi4_3*pi4*rho0*(P0 + rho0) &
!            + 0.5d0*(Nu4 - pi4_3*rho0*Nu2) + 0.5d0*wc2b*(Nu2 - pi4_3*2.d0*rho0))* W0
!        Z(4) = b1+b2+b3
!
!        call LinSys_Cramer(U,Z,a(2,:))

!        ld_r0_reg(1:4) = a(1,1:4) + 0.5d0*a(2,1:4)*r**2
        ld_r0_reg(1:4) = a(1,1:4)
    end function ld_r0_reg

end module
