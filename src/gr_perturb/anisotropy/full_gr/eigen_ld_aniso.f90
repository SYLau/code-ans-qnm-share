module eigen_ld_aniso_mod
    !This module contains several routines to solve for the QNM
    implicit none
    private
    public::pul_sol
    public::scan_Ain_r
    public::get_fmode
    public::brac_eigenmodes_ld
    public::get_eigenmodes_ld
    public::get_ld_var
    type pul_sol
        real(8),allocatable::t(:)
        complex(8),allocatable::x(:,:)
        real(8),allocatable,dimension(:)::p,rho,m,nu,ga,A
        complex(8),allocatable,dimension(:)::H1,K,W,Xp,H0,V,Z

        real(8),allocatable,dimension(:)::rext
        complex(8),allocatable,dimension(:)::Ze, dZe
    end type

    integer :: muller_max_steps = 100
    logical :: print_muller = .false.

    contains
    subroutine scan_Ain_r(l,f1,f2,step,ei_fc,ei_Ain,focus_n,focus_f,focus_range,fi) !Scan the ingoing wave amplitude gamma as a function of frequency (Hz)
        use type_mod,only:pi2
        use puls_ld_aniso_mod,only:integrate_ld
        use io_mod,only:reallocate_a
        use nr_mod,only:sort
        real(8),intent(in)::l,f1,f2
        real(8),intent(in),optional::focus_f(:), focus_range(:)
        integer,intent(in),optional::focus_n
        real(8),intent(in),optional::fi
        integer,intent(in)::step
        complex(8),allocatable,intent(out)::ei_fc(:),ei_Ain(:) !output a list of frequency and ingoing wave amplitude
        real(8)::ei_f(step)
        real(8)::df, f_low, f_up
        real(8)::fimag
        integer::i, j,n,posi,fstep,rstep

        allocate(ei_fc(step),ei_Ain(step))

        fimag = 0.d0
        if (present(fi)) fimag = fi

        rstep=step
        posi=0
        if (present(focus_n) .and. present(focus_f) .and. present(focus_range)) then
            n=size(focus_f)
            if (size(focus_range)/= n) then
                write(*,*)'err: scan_Ain_r focus_range size mismatch'
                stop
            endif
            fstep=step*focus_n/100/n
            if (fstep==0) then
                write(*,*)'err: scan_Ain_r fstep=0'
                stop
            endif
            do i=1,n
                f_low=focus_f(i)-abs(focus_range(i))
                f_up=focus_f(i)+abs(focus_range(i))
                df= (f_up-f_low)/fstep
                forall (j=1:fstep) ei_f(j+posi)=f_low+df*j
                rstep=rstep-fstep
                posi=posi+fstep
            end do
        end if
        if (rstep>0) then
            df=(f2-f1)/rstep
            forall (i=1:rstep) ei_f(i+posi)=f1+df*i
        endif
        call sort(ei_f)
        ei_fc=cmplx(ei_f,fimag,kind=8)
        do i=1,step
            !ei_fc(i+posi)=cmplx(f1+df*i,0.d0,kind=8)
            call integrate_ld(l,pi2*ei_fc(i),.false.,ei_Ain(i))
        end do
    end subroutine

    subroutine get_fmode(l,fguess,fc_save,puls,quit_if_error,ostat) !Search for the position of f-mode along Re(f) and use Muller's method to solve for the QNM
        use type_mod,only:pi2
        use nr_mod,only:mnbrak,golden
        real(8),intent(in)::l,fguess
        logical,intent(in),optional::quit_if_error
        integer,intent(out),optional::ostat
        complex(8),allocatable,intent(out)::fc_save(:)
        type(pul_sol),allocatable,intent(out)::puls(:)
        complex(8)::x1,x2,x3
        real(8)::ax,bx,cx,fa,fb,fc,xm,fm

        allocate(fc_save(1),puls(1))
        ! use golden section search for the minimum point of gamma
        ax=fguess; bx=ax*1.2d0
        call mnbrak(ax,bx,cx,fa,fb,fc,ei_cond_real)
        fm=golden(ax,bx,cx,ei_cond_real,1.d-6,xm)
        x1=cmplx(xm,0.d0,kind=8); x2=x1*1.01d0; x3=x1*0.99d0

        call get_eigenmodes_ld(l,x1,x2,x3,fc_save(1),puls(1),quit_if_error,ostat)

        contains
        function ei_cond_real(x)
            use puls_ld_aniso_mod,only:integrate_ld
            real(8),intent(in) :: x
            real(8) ::ei_cond_real
            complex(8)::yout,w

            w=pi2*cmplx(x,0.d0,kind=8)
            call integrate_ld(l,w,.false.,yout)
            ei_cond_real=abs(yout)
        end function ei_cond_real
    end subroutine get_fmode

    subroutine brac_eigenmodes_ld(l,ei_fc,ei_Ain,fc_save,puls,quit_if_error,ostat)
        !Analyze A_in vs f; use get_eigenmodes_ld to find the roots
        !Locate the minimum and use get_eigenmodes_ld to find the complex root
        !Algorithm not yet optimized
        use type_mod,only:pi2
        use nr_mod,only:mnbrak,golden
        use io_mod,only:reallocate_a
        real(8),intent(in)::l
        complex(8),intent(in)::ei_fc(:),ei_Ain(:)
        logical,intent(in),optional::quit_if_error
        integer,intent(out),optional::ostat
        complex(8),allocatable,intent(out)::fc_save(:)
        type(pul_sol),allocatable,intent(out)::puls(:)
        complex(8)::x1,x2,x3
        real(8),parameter::eps=1.d-1,f_brac=200.d0
        real(8)::dA(size(ei_fc)-1)
        real(8)::ax,bx,cx,fa,fb,fc,xm,fm, A_scale
        integer::i,j,n
        integer,parameter::max_mode=15

        !if (allocated(puls)) deallocate(puls)
        allocate(fc_save(0),puls(0))
        forall (i=1:size(ei_fc)-1) dA(i)=abs(ei_Ain(i+1))-abs(ei_Ain(i))
        do i=1,size(ei_fc)-2
            if (dA(i+1)>0.d0 .and.dA(i)<0.d0) then ! .and. min(abs(ei_Ain(i)),abs(ei_Ain(i+1))) <= eps*sum(abs(ei_Ain))/size(ei_Ain)
                ax=real(ei_fc(i));bx=ei_fc(i+1) !too close
                call mnbrak(ax,bx,cx,fa,fb,fc,ei_cond_real)
                fm=golden(ax,bx,cx,ei_cond_real,1.d-8,xm)

                A_scale=0.d0;n=0 !find a scale for A_in by averaging over a freq within a range of +/- f_brac
                do j=1,size(ei_fc)
                    if (abs(ei_fc(j)-xm)<f_brac) then
                        A_scale=A_scale+abs(ei_Ain(j))
                        n=n+1
                    end if
                end do
                if (n/=0) A_scale=A_scale/n
                print*, 'A_scale = ', A_scale

                if (abs(ei_cond_real(xm))>eps*A_scale) cycle !reject fake minima
                if (size(fc_save)>max_mode) then
                    write(*,*)'warning: get_eigenmodes_ld exceed max. no. of modes: ',max_mode
                    exit
                endif

                fc_save=reallocate_a(fc_save,size(fc_save)+1)
                puls=reallocate_pul_sol_a(puls,size(puls)+1)
                !x1=ei_fc(i); x2=ei_fc(i+1); x3=(ei_fc(i)+ei_fc(i+1))/2.d0
                x1=cmplx(xm,0.d0,kind=8); x2=x1*1.001d0; x3=x1*0.999d0
                call get_eigenmodes_ld(l,x1,x2,x3,fc_save(size(fc_save)),puls(size(puls)),quit_if_error,ostat)
            endif
        end do
        contains
        function ei_cond_real(x)
            use puls_ld_aniso_mod,only:integrate_ld
            real(8),intent(in) :: x
            real(8) ::ei_cond_real
            complex(8)::yout,w
            w=pi2*cmplx(x,0.d0,kind=8)
            call integrate_ld(l,w,.false.,yout)
            ei_cond_real=abs(yout)
        end function ei_cond_real
    end subroutine brac_eigenmodes_ld

    subroutine get_eigenmodes_ld(l,x1,x2,x3,fc_save,puls_save,quit_if_error,ostat) !Use Muller's method to solve for complex roots
        use type_mod,only:pi2
        use puls_ld_aniso_mod,only:integrate_ld,t_pu,x_pu
        use puls_ld_aniso_mod,only:ext
        use muller_mod,only:rtmuller
        real(8),intent(in)::l
        complex(8),intent(in)::x1,x2,x3
        logical,intent(in),optional::quit_if_error
        integer,intent(out),optional::ostat
        complex(8),intent(out)::fc_save
        type(pul_sol),intent(out)::puls_save
        complex(8)::dump
        real(8)::eps=1.d-11
            fc_save=rtmuller(ei_cond, x1, x2, x3, eps, quit_if_error,ostat,maxit=muller_max_steps) !use Muller's method to locate the mode frequency
            if (present(ostat)) then
                if (ostat/=0) return
            end if
            call integrate_ld(l,pi2*fc_save,.true.,dump) !write solution at the mode frequency
            !puls=reallocate_pul_sol_a(puls,size(puls)+1)
            puls_save%t=t_pu
            puls_save%x=x_pu
            deallocate(t_pu,x_pu)

            puls_save%rext=ext%r
            puls_save%Ze=ext%Ze
            puls_save%dZe=ext%dZe
            deallocate(ext%r,ext%Ze, ext%dZe)

        contains
        function ei_cond(x)
            complex(8),intent(in) :: x
            complex(8)::ei_cond
            complex(8)::yout,w
            w=pi2*x
            call integrate_ld(l,w,.false.,yout)
            ei_cond=yout
            if (print_muller) print*,w/pi2,abs(yout)
        end function ei_cond
    end subroutine get_eigenmodes_ld

    subroutine get_ld_var(l,fc,puls) !Convert the dependent variables in the ODEs into the corresponding variables in LD formalism for output
        use type_mod,only:G,c,pi2
        use bg_mod,only:get_bg_r,get_tov_metric
        use ld_eq_aniso_mod,only:get_H0_V_Z,perturb_par
        real(8),intent(in)::l
        complex(8),intent(in)::fc(:)
        type(perturb_par)::par
        type(pul_sol),allocatable,intent(inout)::puls(:)
        real(8)::gc2,gc4,p,rho,m,nu,ga,A,eLam,dLam,dnu,eNu
        complex(8)::H0,V,Z,wc2
        integer::i,j,isave, size_t
        !if (.not. allocated(puls)) pause 'err: get_ld_var puls not allocated'
        gc2=G/c**2; gc4=G/c**4
        if (size(fc)/=size(puls)) error stop 'err: get_ld_var size fc and puls mismatch'

        do i=1,size(puls)
            if (.not. allocated(puls(i)%t)) error stop 'err: get_ld_var puls(i)%t not allocated'

            size_t=size(puls(i)%t)
            wc2=(pi2*fc(i)/c)**2
            isave=1
            if (allocated(puls(i)%p)) deallocate(puls(i)%p)
            if (allocated(puls(i)%rho)) deallocate(puls(i)%rho)
            if (allocated(puls(i)%m)) deallocate(puls(i)%m)
            if (allocated(puls(i)%nu)) deallocate(puls(i)%nu)
            if (allocated(puls(i)%ga)) deallocate(puls(i)%ga)
            if (allocated(puls(i)%A)) deallocate(puls(i)%A)

            if (allocated(puls(i)%H1)) deallocate(puls(i)%H1)
            if (allocated(puls(i)%K)) deallocate(puls(i)%K)
            if (allocated(puls(i)%W)) deallocate(puls(i)%W)
            if (allocated(puls(i)%Xp)) deallocate(puls(i)%Xp)
            if (allocated(puls(i)%H0)) deallocate(puls(i)%H0)
            if (allocated(puls(i)%V)) deallocate(puls(i)%V)
            if (allocated(puls(i)%Z)) deallocate(puls(i)%Z)
            allocate(puls(i)%p, puls(i)%rho, puls(i)%m, puls(i)%nu, puls(i)%ga, puls(i)%A &
            ,source = puls(i)%t)
            allocate(puls(i)%H1(size_t),puls(i)%K(size_t),puls(i)%W(size_t),puls(i)%Xp(size_t) &
            ,puls(i)%H0(size_t),puls(i)%V(size_t),puls(i)%Z(size_t))

            puls(i)%H1(:)=puls(i)%x(1,:)
            puls(i)%K(:)=puls(i)%x(2,:)
            puls(i)%W(:)=puls(i)%x(3,:)
            puls(i)%Xp(:)=puls(i)%x(4,:)

            do j=1,size_t
                call get_bg_r(puls(i)%t(j),p,rho,m,nu,ga,A,isave)
                call get_tov_metric(puls(i)%t(j),p,rho,m,nu,eLam,dLam,dnu,eNu)

                puls(i)%p(j) = p
                puls(i)%rho(j) = rho
                puls(i)%m(j) = m
                puls(i)%nu(j) = nu
                puls(i)%ga(j) = ga
                puls(i)%A(j) = A

!                p=gc4*p; rho=gc2*rho; m=gc2*m !unit conversion

                par%l=l; par%wc2=wc2; par%p=gc4*p; par%rho=gc2*rho; par%m=gc2*m; par%nu=nu; par%ga=ga
                par%eLam=eLam; par%dLam=dLam; par%dnu=dnu; par%eNu=eNu; par%A=A

                call get_H0_V_Z(par,puls(i)%t(j),puls(i)%x(:,j),H0,V,Z)

                puls(i)%H0(j)=H0
                puls(i)%V(j)=V
                puls(i)%Z(j)=Z

            end do

        end do

    end subroutine get_ld_var


!==================================================================================
!    local utility
    function reallocate_pul_sol_a(p,n)
        type(pul_sol), dimension(:), allocatable :: p, reallocate_pul_sol_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_pul_sol_a(n),stat=ierr)
        if (ierr /= 0) then
            write(*,*)'reallocate_pul_sol_a: problem in attempt to allocate memory'
            stop
        endif
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_pul_sol_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_pul_sol_a
end module

