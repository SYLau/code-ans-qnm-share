module tov_aniso_pt_mod
    !This module solves the TOV equations when provided the EOS in the format specified by the interface below [see the subroutine solve_tov_aniso_pt]
    implicit none
    private
    public::set_eps
    public::tov_aniso_pt_bg
    public::solve_tov_aniso_pt
    public::Sch
    public::Sch_barotropic, Sch_general

    type tov_aniso_pt_bg
        real(8),allocatable,dimension(:)::r,p,rho,m,nu,ga,A,s !Schwarzchild discriminant "A" default to be zero (barotrope)
        integer,allocatable,dimension(:)::pt_pos
        integer::n
    end type

    real(8)::set_eps=1.d-12
    real(8),allocatable,private:: t_bg(:),x_bg(:,:)
    !type(tov_aniso_pt_bg),allocatable,protected,public:: tovs(:)
    !integer,allocatable,protected,public:: pt_i(:)

    real(8),parameter,private::eps_r0=1.d-7 !7.15d-7 !initial radius

    interface
        function Sch_func(r,p,rho,m,ga) result(A)
            real(8),dimension(:),intent(in)::r,p,rho,m,ga
            real(8),dimension(size(r))::A
        end function Sch_func
    end interface

    procedure(Sch_func),pointer:: Sch => Sch_barotropic


contains

    subroutine solve_tov_aniso_pt(Pc,pt,eps_Pm,eos_rho,eos_ga_r,sigma,tovs) !solve TOV eq for models with multiple 1st order phase transitions; store results in tovs (derived data type)
        use type_mod,only:G,pi,c
        use nr_mod,only:polint
        use io_mod,only:reallocate_a
        implicit none
        real(8),intent(in)::Pc,pt(:),eps_Pm
        type(tov_aniso_pt_bg),intent(out)::tovs
        real(8)::r0,P0,m0,nu0,rc,rhoc,R2,Pm,h0,nuR
        real(8)::err
        real(8),parameter::eps_small=5.d-16
        integer::i,j,nold
        interface
            function eos_rho(p) !input p return rho
                real(8),intent(in)::p
                real(8)::eos_rho
            end function eos_rho
            function eos_ga_r(p)
                real(8),intent(in)::p
                real(8)::eos_ga_r
            end function eos_ga_r
            function sigma(p,rho,m,r)
                real(8),intent(in)::p,rho,m,r
                real(8)::sigma
            end function sigma
        end interface

        rhoc= eos_rho(Pc)
        R2=sqrt(3.d0*Pc/(2.d0*pi*G*(rhoc+Pc/c**2)*(rhoc+3.d0*Pc/c**2))) !Taylor's expansion of dP/dr; see e.g. Kruger et al. 2015
        rc= eps_r0*R2

        allocate(tovs%r(0))
        allocate(tovs%p,tovs%rho,tovs%m,tovs%nu,tovs%ga,tovs%A,tovs%s,source=tovs%r)
        allocate(tovs%pt_pos(size(pt)))
        i=1
        r0=rc
        P0=Pc*(1.d0-(rc/R2)**2)
        m0=4.d0*pi/3.d0*rhoc*rc**3
        nu0=2.d0*pi/3.d0*G/c**2*(rhoc+3.d0*Pc/c**2)*rc**2
        h0=rc

        do
            if (size(pt)-i+1>0) Pm=pt(i)
            if (size(pt)-i+1==0) Pm=eps_Pm*Pc

            call integrate_tov_aniso(P0,m0,nu0,r0,Pm,h0,eos_rho,sigma)

            nold=size(tovs%r)
            tovs%r=reallocate_a(tovs%r,nold+size(t_bg))
            tovs%p=reallocate_a(tovs%p,nold+size(t_bg))
            tovs%rho=reallocate_a(tovs%rho,nold+size(t_bg))
            tovs%m=reallocate_a(tovs%m,nold+size(t_bg))
            tovs%nu=reallocate_a(tovs%nu,nold+size(t_bg))
            tovs%ga=reallocate_a(tovs%ga,nold+size(t_bg))
            tovs%A=reallocate_a(tovs%A,nold+size(t_bg))
            tovs%s=reallocate_a(tovs%s,nold+size(t_bg))
            if (size(pt)-i+1>0) tovs%pt_pos(i)=nold+size(t_bg)

!            do j=1,size(t_bg)-1
            do j=1,size(t_bg)
                tovs%r(j+nold)=t_bg(j)
                tovs%p(j+nold)=x_bg(1,j)
                tovs%rho(j+nold)=x_bg(2,j)
                tovs%m(j+nold)=x_bg(3,j)
                tovs%nu(j+nold)=x_bg(4,j)
                tovs%ga(j+nold)=eos_ga_r(x_bg(1,j))
!                tovs%A(j+nold)=0.d0
                tovs%s(j+nold)=sigma(G/c**4*tovs%p(j+nold), G/c**2*tovs%rho(j+nold), G/c**2*tovs%m(j+nold), tovs%r(j+nold))
            end do
            tovs%n=size(tovs%r)

            !interpolation for last data point
!            tovs%n=size(tovs%r)
!            tovs%P(tovs%n)=Pm
!            tovs%rho(tovs%n)=eos_rho(Pm*(1.d0+eps_Pm))
!            tovs%ga(tovs%n)=eos_ga_r(Pm*(1.d0+eps_Pm))
!
!            call polint(x_bg(1,size(t_bg)-1:size(t_bg)),t_bg(size(t_bg)-1:size(t_bg)),Pm,tovs%r(tovs%n),err)
!            if (tovs%r(tovs%n)==tovs%r(tovs%n-1)) then
!                tovs%r(tovs%n)=tovs%r(tovs%n)*(1+eps_small) ! prevent underflow of independent variable r
!                if (tovs%r(tovs%n)==tovs%r(tovs%n-1)) then
!                    print*, 'err: solve_tov_aniso_pt r underflow, increase eps_small'
!                    stop
!                end if
!            endif
!            call polint(x_bg(1,size(t_bg)-1:size(t_bg)),x_bg(3,size(t_bg)-1:size(t_bg)),Pm,tovs%m(tovs%n),err)
!            call polint(x_bg(1,size(t_bg)-1:size(t_bg)),x_bg(4,size(t_bg)-1:size(t_bg)),Pm,tovs%nu(tovs%n),err)
!            tovs%s(tovs%n)=sigma( G/c**4*tovs%p(tovs%n), G/c**2*tovs%rho(tovs%n), G/c**2*tovs%m(tovs%n), tovs%r(tovs%n) )

            tovs%A(nold+1:size(tovs%r)) = Sch(tovs%r(nold+1:size(tovs%r)),tovs%p(nold+1:size(tovs%r)) &
            ,tovs%rho(nold+1:size(tovs%r)),tovs%m(nold+1:size(tovs%r)),tovs%ga(nold+1:size(tovs%r)))

            deallocate(t_bg,x_bg)
            if (size(pt)-i+1==0) exit
            i=i+1
            r0=tovs%r(tovs%n)
            P0=tovs%P(tovs%n)*(1.d0-eps_Pm)
            m0=tovs%m(tovs%n)
            nu0=tovs%nu(tovs%n)
        enddo

        nuR=0.5d0*log(1.d0-2.d0*G*tovs%m(tovs%n)/tovs%r(tovs%n)/c**2)-tovs%nu(tovs%n)
        tovs%nu=tovs%nu+nuR
!        tovs%A=0 !Barotrope
    end subroutine

    !============================================================================
    ! Use RK4 to integrate TOV equation
    ! Solution is stored in module variable t_bg, x_bg
    ! Integration ends when P <= Pm
    !============================================================================
    subroutine integrate_tov_aniso(P0,m0,nu0,r0,Pm,h0,eos_rho,sigma)
        use type_mod,only:G,pi,c
!        use odeSolver_tov
        use odeSolver_mod, only: rk45ad
        use odeSolver_mod, only: tp, xp
        implicit none
        real(8),intent(in)::P0,m0,nu0,r0,Pm
        real(8),intent(inout)::h0
        real(8)::t,x(1:3)
        real(8)::h,emax
        integer:: itmax=100, iflag
        integer::i
        interface
            function eos_rho(p)
                real(8),intent(in)::p
                real(8)::eos_rho
            end function eos_rho
            function sigma(p,rho,m,r)
                real(8),intent(in)::p,rho,m,r
                real(8)::sigma
            end function sigma
        end interface

        emax=set_eps

        t=r0
        x(1)=P0
        x(2)=m0
        x(3)=nu0
        h=h0

!        call rk45ad_ob(tov_eq,tov_cond,t,x,h,itmax,emax,iflag)
        call rk45ad(tov_eq,t,x,h,r0*1.e10,emax,tcond=tov_cutoff,openbd=.true.)

        if (allocated(t_bg)) deallocate(t_bg)
        if (allocated(x_bg)) deallocate(x_bg)
        allocate(t_bg(size(tp)),x_bg(4,size(tp)))
        t_bg=tp
        x_bg(1,:)=xp(1,:)
        x_bg(2,1) = eos_rho(xp(1,1)*(1-epsilon(eos_rho(xp(1,1)))))
        do i=2,size(tp)-1
            x_bg(2,i)=eos_rho(xp(1,i))
        enddo
        x_bg(2,size(tp)) = eos_rho(xp(1,size(tp))*(1+epsilon(eos_rho(xp(1,size(tp))))))
        x_bg(3:4,:)=xp(2:3,:)
        deallocate(tp,xp)

        h0=h

    contains
        subroutine tov_eq(t,x,f)
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            real(8)::r,rho,p,m,nu
            r=t;p=x(1);m=x(2);nu=x(3);rho=eos_rho(p)
            f(1)=-G/r**2*(rho+p/c**2)*(m+4.d0*pi*r**3*P/c**2)/(1.d0-2.d0*G*m/c**2/r) &
            -2*sigma(G/c**4*p,G/c**2*rho,G/c**2*m,r)/r/(G/c**4)
            f(2)=4.d0*pi*rho*r**2
            f(3)=G/c**2/r**2*(m+4.d0*pi*r**3*P/c**2)/(1.d0-2.d0*G*m/c**2/r)
        endsubroutine tov_eq
!        function tov_cond(t,x)
!            real(8), intent(in):: t, x(:)
!            logical:: tov_cond
!            tov_cond = (x(1)<=Pm)
!        endfunction tov_cond
        function tov_cutoff(t,x) result(res)
            real(8), intent(in):: t, x(:)
            real(8):: res
            res = x(1)-Pm
        endfunction tov_cutoff
    end subroutine integrate_tov_aniso

    !============================================================================
    ! Schwarzschild discriminant
    ! Barotropic case
    !============================================================================
    function Sch_barotropic(r,p,rho,m,ga) result(A)
        real(8),dimension(:),intent(in)::r,p,rho,m,ga
        real(8),dimension(size(r))::A
        A=0
    end function Sch_barotropic

    !============================================================================
    ! Schwarzschild discriminant
    ! General case with EOS table
    !============================================================================
    function Sch_general(r,p,rho,m,ga) result(A)
        use type_mod,only: G,pi,c
        real(8),dimension(:),intent(in)::r,p,rho,m,ga
        real(8),dimension(size(r))::A
        real(8),dimension(size(r))::dpr,ga0
        integer::n,i
        n=size(r)
        dpr = -G/r**2*(rho+p/c**2)*(m+4*pi*r**3*P/c**2)/(1-2*G*m/c**2/r)
        do i=1,n-1
            ga0(i)=(log(p(i+1))-log(p(i)))/(log(rho(i+1))-log(rho(i)))*(1+p(i)/rho(i)/c**2)
        end do
        ga0(n)=(log(p(n))-log(p(n-1)))/(log(rho(n))-log(rho(n-1)))*(1+p(n)/rho(n)/c**2)
        A = (1/ga0-1/ga)/p*dpr
    end function Sch_general

end module

