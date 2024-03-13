module puls_ld_aniso_mod
    !This module contains the integration of the 5 independent solutions of LD formalism
    !as well as the matching of them at the matching position
    !For models without 1st order phase transition, the matching point is half the radius of the star
    !For models with PT, it matches at half of the width of the region closest to the surface
    use type_mod,only:wp
    implicit none
    private
    public::integrate_ld
    public::ode_eps,ode_eps_hmin
    public::ug_steps
    public::n_initial

    public::exterior_method

    public::t_pu,x_pu
    public::ext
    public::psave

    !Integration parameters
    real(wp)::ode_eps=1.e-9_wp
    real(wp)::ode_eps_hmin=0._wp
!    real(wp)::ode_eps_hmin=1.d-13
    integer::ug_steps=3000
    integer::n_initial=1

    !Exterior problem: Zerilli's method r_infinity
    real(wp)::rinf_Zerilli_default = 25._wp
    !Exterior problem method option (1: Zerilli; 2. Leaver)
    integer::exterior_method = 1

    !Note (07/08/2022):
    !ode_eps_hmin is needed to set a minimum hmin in rk45ad, or else the integral near surface never converge (probably at low frequency)
    !also when mode frequency has negative imaginary part, it would cause problem
    !need to use uniform grid to get i-mode; mode function near surface for i-mode looks discrete
    !the old code has smoothened the gamma by taking numerical derivative of the 4th order interpolated P and rho

    real(8),allocatable::t_pu(:)
    complex(8),allocatable::x_pu(:,:)

    type ext_sol
        real(wp),allocatable,dimension(:)::r
        complex(wp),allocatable,dimension(:)::Ze, dZe
    end type ext_sol

    type (ext_sol)::ext


    type pul_sol
        real(8),allocatable::t(:)
        complex(8),allocatable::x(:,:)
    end type pul_sol

    type(pul_sol),dimension(5)::psave

contains
    subroutine integrate_ld(l,w,savedata,yout ,method) ! Input ell, omega (s^{-1}), returns the ingoing wave amplitude gamma as yout
        use odeSolver_mod,only:rk45ad,tp,xp
        use bg_mod,only:tovs=>bgs,pt_i
        use io_mod,only:reallocate_a
        implicit none
        complex(8),intent(in)::w
        real(8),intent(in)::l
        logical,intent(in)::savedata
        complex(8),intent(out)::yout

        character(len=*),optional,intent(in)::method

        character(len=:),allocatable::ode_method
        real(8)::t,x(8),tb,xv(8,5)
        complex(8)::xvc(4,5),H0,K,Wm
        integer::isave,i,nold,isol
        real(8)::emax
        complex(8)::coef(5)

        ode_method='lsoda'
        if (present(method)) ode_method=method

        emax=ode_eps
!        eps_hmin=ode_eps_hmin
        loop_frd: do isol=1,2
            isave=1 !save r position for faster interpolation with get_bg
            t=tovs(n_initial)%r
            x=bc_initial(isol,(1.d0,0.d0))
            i=1
            !save partial solutions
            if (allocated(psave(isol)%t)) deallocate(psave(isol)%t)
            if (allocated(psave(isol)%x)) deallocate(psave(isol)%x)
            allocate(psave(isol)%t(0),psave(isol)%x(4,0))
            do
                if (size(pt_i)-i+1>0) tb=tovs(pt_i(i))%r
                if (size(pt_i)-i+1==0) tb=(tovs(size(tovs))%r+t)/2.d0

                call integrate_ode(ode_method,t,x,tb,emax)

                !save partial solutions for testing
                nold=size(psave(isol)%t)
                psave(isol)%t=reallocate_a(psave(isol)%t,nold+size(tp))
                psave(isol)%x=reallocate_a(psave(isol)%x,4,nold+size(tp))
                psave(isol)%t(nold+1:size(psave(isol)%t))=tp(1:size(tp))
                psave(isol)%x(1:4,nold+1:size(psave(isol)%t))=cmplx(xp(1:4,1:size(tp)),xp(5:8,1:size(tp)),kind=wp)

                if (size(pt_i)-i+1==0) exit
                t=tovs(pt_i(i)+1)%r
                x=x !LD variables are continuous across fluid-fluid interface
                i=i+1
            enddo
            xv(:,isol)=x

        enddo loop_frd

        !back integration
        loop_brd: do isol=3,5
            isave=size(tovs)/2
            t=tovs(size(tovs))%r
            x=bc_surface(isol,(1.d0,0.d0))
            !save partial solutions
            if (allocated(psave(isol)%t)) deallocate(psave(isol)%t)
            if (allocated(psave(isol)%x)) deallocate(psave(isol)%x)
            allocate(psave(isol)%t(0),psave(isol)%x(4,0))

            call integrate_ode(ode_method,t,x,tb,emax)

            !save partial solutions for testing
            nold=size(psave(isol)%t)
            psave(isol)%t=reallocate_a(psave(isol)%t,nold+size(tp))
            psave(isol)%x=reallocate_a(psave(isol)%x,4,nold+size(tp))
            psave(isol)%t(nold+1:size(psave(isol)%t))=tp(size(tp):1:-1)!reverse order
            psave(isol)%x(1:4,nold+1:size(psave(isol)%t))=cmplx(xp(1:4,size(tp):1:-1),xp(5:8,size(tp):1:-1),kind=wp)
            xv(:,isol)=x
        enddo loop_brd

        !solve for surface sol
        xvc(1:4,1:5)=cmplx(xv(1:4,:),xv(5:8,:),kind=wp)
        coef=get_sol_coef(xvc)
        !solve Zerilli eq in Vacuum
        t=tovs(size(tovs))%r

        x=bc_surface(3,coef(3))+bc_surface(4,coef(4))+bc_surface(5,coef(5))
        Wm=dot_product(conjg(xvc(3,3:5)),coef(3:5)) !necessary?
        H0=get_H0(t,x)/Wm
        K=coef(4)/Wm

        x(1:4)=bc_vac_initial(H0,K)

        solve_exterior: block
            use type_mod,only:G,c
            real(8)::h,hmin
            real(8)::eps_hmin

            eps_hmin=ode_eps_hmin
            if (exterior_method == 1) then
                ! Zerilli's method
                tb = rinf_Zerilli_default*c/real(w)
                h=(tb-t)/3000
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_vac,t,x(1:4),h,tb,emax,savedata,hmin=hmin)
                yout=bc_infty(t,x(1:4)) !obtain the ingoing Zerilli function amplitude as the output

                if (savedata) then
                    ext%r = tp
                    ext%Ze = cmplx(xp(1,:), xp(3,:),kind=wp)*Wm !Remove normalization so that the saved solution has the same normalization (unnormalized) as the interior saved solutions
                    ext%dZe = cmplx(xp(2,:), xp(4,:),kind=wp)*Wm
                end if
            elseif (exterior_method == 2) then
                ! Leaver's method
                if (4*G*tovs(size(tovs))%m/c**2 > tovs(size(tovs))%r) then
                    tb=4*G*tovs(size(tovs))%m/c**2
                    h=(tb-t)/3000
                    hmin=(tb-t)*eps_hmin
                    call rk45ad(ode_vac,t,x(1:4),h,tb,emax,savedata,hmin=hmin)
                else
                    tb=tovs(size(tovs))%r
                    t=tb
                end if
                yout=bc_leaver(t,x(1:4))

                if (savedata) then
                    ext%r = tp
                    ext%Ze = cmplx(xp(1,:), xp(3,:),kind=wp)*Wm
                    ext%dZe = cmplx(xp(2,:), xp(4,:),kind=wp)*Wm
                end if
            else
                error stop 'err: integrate_ld exterior_method not specified'
            end if
        end block solve_exterior

        !integrate once more with correct BCs and save solution
        if (savedata) then
            if (allocated(t_pu)) deallocate(t_pu)
            if (allocated(x_pu)) deallocate(x_pu)
            allocate(t_pu(0),x_pu(4,0))
            isave=1 !save r position for faster interpolation with get_bg
            t=tovs(n_initial)%r
            x=bc_initial(1,coef(1))+bc_initial(2,coef(2)) !The correct BC
            i=1
            do
                if (size(pt_i)-i+1>0) tb=tovs(pt_i(i))%r
                if (size(pt_i)-i+1==0) tb=(tovs(size(tovs))%r+t)/2.d0

                call integrate_ode(ode_method,t,x,tb,emax)

                !save solutions
                nold=size(t_pu)
                t_pu=reallocate_a(t_pu,nold+size(tp))
                x_pu=reallocate_a(x_pu,4,nold+size(tp))
                t_pu(nold+1:size(t_pu))=tp(1:size(tp))
                x_pu(1:4,nold+1:size(t_pu))=cmplx(xp(1:4,1:size(tp)),xp(5:8,1:size(tp)),kind=wp)

                if (size(pt_i)-i+1==0) exit
                t=tovs(pt_i(i)+1)%r
                x=x !LD variables are continuous
                i=i+1
            enddo
            isave=size(tovs)/2
            t=tovs(size(tovs))%r
            x=bc_surface(3,coef(3))+bc_surface(4,coef(4))+bc_surface(5,coef(5)) !The correct BC

            call integrate_ode(ode_method,t,x,tb,emax)

            nold=size(t_pu)
            t_pu=reallocate_a(t_pu,nold+size(tp))
            x_pu=reallocate_a(x_pu,4,nold+size(tp))
            t_pu(nold+1:size(t_pu))=tp(size(tp):1:-1)!reverse order
            x_pu(1:4,nold+1:size(t_pu))=cmplx(xp(1:4,size(tp):1:-1),xp(5:8,size(tp):1:-1),kind=wp)
        endif
    contains

        !Interface connecting with ODEs and BCs
        subroutine integrate_ode(ode_method,t,x,tb,emax)
            use odeSolver_mod, only:rk4ug
            use odepdr_mod,only:odelsoda
            character(len=*),intent(in)::ode_method
            real(wp),intent(inout)::t,x(:)
            real(wp),intent(in)::tb,emax
            real(wp)::h,hmin
            real(wp)::eps_hmin

            eps_hmin=ode_eps_hmin
            if (ode_method=='rk45') then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_stellar,t,x,h,tb,emax,.true.,hmin=hmin)
            else if (ode_method=='uniform') then
                call rk4ug(ode_stellar,t,x,ug_steps,tb,.true.)
            else if (ode_method=='lsoda') then
                call odelsoda(ode_stellar,t,x,tb,emax/2.e2,t_len=3000)
            else
                error stop 'err: ode_method not valid'
            endif
        end subroutine integrate_ode


        function bc_initial(isol,coef)
            use type_mod,only:G,c
            use ld_eq_aniso_mod,only:ld_r0_reg
            complex(8),intent(in)::coef
            integer,intent(in)::isol
            real(8)::P0,rho0,nu0,ga0,r0,bc_initial(8)
            P0 = G/c**4*tovs(n_initial)%p
            rho0 = G/c**2*tovs(n_initial)%rho
            Nu0 = tovs(n_initial)%nu
            ga0 = tovs(n_initial)%ga
            r0 = tovs(n_initial)%r
            bc_initial(1:4)=real(coef * ld_r0_reg(isol,P0,rho0,nu0,ga0,r0,l,w/c))
            bc_initial(5:8)=imag(coef * ld_r0_reg(isol,P0,rho0,nu0,ga0,r0,l,w/c))
        end function

        function bc_surface(isol,coef)
            complex(8),intent(in)::coef
            integer,intent(in)::isol
            real(8)::bc_surface(8)
            bc_surface=0.d0
            if (isol==3) bc_surface(1)=real(coef) !3 independent surface sol.
            if (isol==3) bc_surface(5)=imag(coef)
            if (isol==4) bc_surface(2)=real(coef)
            if (isol==4) bc_surface(6)=imag(coef)
            if (isol==5) bc_surface(3)=real(coef)
            if (isol==5) bc_surface(7)=imag(coef)
        end function

        function get_sol_coef(xc)
            use linsys_mod,only:LinSys_Cramer
            complex(8),intent(in)::xc(:,:)
            complex(8)::get_sol_coef(5)
            complex(8)::U(4,4),coefc(4)
            U(:,1)=-xc(:,2); U(:,2:4)=xc(:,3:5)
            call LinSys_Cramer(U,xc(:,1),coefc)
            get_sol_coef(1)=(1.d0,0.d0)
            get_sol_coef(2:5)=coefc(1:4)
        end function

        subroutine ode_stellar(t,x,f)
            use type_mod,only:G,c
            use bg_mod,only:get_tov_metric,get_bg_r
            use ld_eq_aniso_mod,only:ld_eq, perturb_par
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            type(perturb_par)::par
            real(8)::gc2,gc4,p,rho,m,nu,ga,A,eLam,dLam,dnu,eNu
            complex(8)::xc(4),fc(4),wc2

            gc2=G/c**2
            gc4=G/c**4
            wc2=(w/c)**2
            xc(1:4) = cmplx(x(1:4),x(5:8),kind=8)

            call get_bg_r(t,p,rho,m,nu,ga,A,isave)
            call get_tov_metric(t,p,rho,m,nu,eLam,dLam,dnu,eNu)

            par%l=l; par%wc2=wc2; par%p=gc4*p; par%rho=gc2*rho; par%m=gc2*m; par%nu=nu; par%ga=ga
            par%A=gc2*A
            par%eLam=eLam; par%dLam=dLam; par%dnu=dnu; par%eNu=eNu
            call ld_eq(par,t,xc,fc)

            f(1:4)=real(fc(1:4))
            f(5:8)=imag(fc(1:4))
        end subroutine

        function get_H0(t,x)
            use type_mod,only:G,c
            use bg_mod,only:get_tov_metric,get_bg_r
            use ld_eq_aniso_mod,only:get_H0_V_Z,perturb_par
            real(8),intent(in)::t,x(8)
            complex(8)::get_H0,H0,V,Z
            type(perturb_par)::par
            real(8)::gc2,gc4,p,rho,m,nu,ga,A,eLam,dLam,dnu,eNu
            complex(8)::xc(4),wc2

            gc2=G/c**2
            gc4=G/c**4
            wc2=(w/c)**2
            xc(1:4) = cmplx(x(1:4),x(5:8),kind=8)
            call get_bg_r(t,p,rho,m,nu,ga,A,isave)
            call get_tov_metric(t,p,rho,m,nu,eLam,dLam,dnu,eNu)

            par%l=l; par%wc2=wc2; par%p=gc4*p; par%rho=gc2*rho; par%m=gc2*m; par%nu=nu; par%ga=ga
            par%A=gc2*A
            par%eLam=eLam; par%dLam=dLam; par%dnu=dnu; par%eNu=eNu

            call get_H0_V_Z(par,t,xc,H0,V,Z)
            get_H0=H0
        end function

        function bc_vac_initial(H0,K)
            use type_mod,only:G,c
            use zerilli_mod,only:zerilli_initial
            complex(8),intent(in)::H0,K
            real(8)::m,r,bc_vac_initial(4)
            m = G/c**2*tovs(size(tovs))%m
            r = tovs(size(tovs))%r
            bc_vac_initial(1:2)=real(zerilli_initial(m,r,H0,K,l,w))
            bc_vac_initial(3:4)=imag(zerilli_initial(m,r,H0,K,l,w))
        end function

        subroutine ode_vac(t,x,f)
            use type_mod,only:G,c
            use zerilli_mod,only:zerilli_eq
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            complex(8)::wc2,xc(2),fc(2)
            real(8)::gc2,m0

            gc2=G/c**2
            wc2=(w/c)**2
            m0=gc2*tovs(size(tovs))%m
            xc(1:2) = cmplx(x(1:2),x(3:4),kind=8)
            call zerilli_eq(l,wc2,m0,t,xc,fc)
            f(1:2)=real(fc(1:2))
            f(3:4)=imag(fc(1:2))
        end subroutine

        function bc_infty(t,x)
            use type_mod,only:G,c
            use zerilli_mod,only:zerilli_asym_coef
            real(8),intent(in)::t,x(:)
            real(8)::m0
            complex(8)::xc(2),wc,be,ga,bc_infty
            xc(1:2)=cmplx(x(1:2),x(3:4),kind=8)
            m0=G/c**2*tovs(size(tovs))%m
            wc=w/c
            call zerilli_asym_coef(l,m0,t,wc,xc,be,ga)
            bc_infty=ga
        end function

        function bc_leaver(t,x) result(res)
            use type_mod,only:G,c
            use leaver_mod,only:RW_Zerilli
            use leaver_mod,only:leaver_solve
            real(wp),intent(in)::t,x(:)
            complex(wp)::res
            complex(wp)::Z(2),Q(2)
            complex(wp)::wc
            real(wp)::m,r
            r=t
            m=G/c**2*tovs(size(tovs))%m
            wc=w/c
            Z=cmplx(x(1:2),x(3:4),kind=wp)
            Q=RW_Zerilli(l,r,m,wc,Z)
            res = leaver_solve(l,r,m,wc,Q)
        end function

    end subroutine integrate_ld

end module
