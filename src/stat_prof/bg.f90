module bg_mod
    !This module contains routines to load in the background solution
    !and routines to interpolate the solutions when given the position r
    !Default interpolation scheme: 4th order polynomial interpolation
    !Also contains the interpolation along the grid x (x=Ln(r/P) ) [ref McDermott 1988]
    implicit none
    private
    public::load_bgs,bg,bgs,pt_i,gx
    public::get_tov_metric,get_tov_metric2,get_bg_r,get_bg_x,get_bg_newt_r,get_x
    public::load_bgs_mu,get_bg_solid_r

    interface load_bgs
        module procedure loadbgs_tov, loadbgs_newt
    end interface load_bgs

    type bg
        real(8)::r,p,rho,m,nu,ga,A
        real(8)::mu,chic !shear modulus; d mu /d P
    end type

    type(bg),allocatable,protected:: bgs(:)
    integer,allocatable,protected:: pt_i(:)

    real(8),allocatable::gx(:)

    real(8),parameter::torr_ri=0.8d0
    real(8),parameter::torr_rf=1.d-15
contains

    subroutine loadbgs_tov(r,p,rho,m,nu,ga,A,pt_index)
        real(8),dimension(:),intent(in)::r,p,rho,m,nu,ga,A
        integer,intent(in)::pt_index(:)
        integer::i
        if (allocated(bgs)) deallocate(bgs)
        if (allocated(pt_i)) deallocate(pt_i)
        allocate(bgs(size(r)),pt_i(size(pt_index)))
        forall (i=1:size(bgs)) bgs(i)=bg(r(i),p(i),rho(i),m(i),nu(i),ga(i),A(i),0.d0,0.d0)
        pt_i=pt_index
    end subroutine

    subroutine loadbgs_newt(r,p,rho,m,ga,A,pt_index)
        real(8),dimension(:),intent(in)::r,p,rho,m,ga,A
        integer,intent(in)::pt_index(:)
        integer::i
        if (allocated(bgs)) deallocate(bgs)
        if (allocated(pt_i)) deallocate(pt_i)
        allocate(bgs(size(r)),pt_i(size(pt_index)))
        forall (i=1:size(bgs)) bgs(i)=bg(r(i),p(i),rho(i),m(i),0.d0,ga(i),A(i),0.d0,0.d0)
        pt_i=pt_index
    end subroutine

    subroutine load_bgs_mu(mu,chic)
        real(8),dimension(:),intent(in)::mu,chic

        bgs%mu=mu
        bgs%chic=chic
    end subroutine load_bgs_mu

    subroutine get_tov_metric(r,p,rho,m,nu,eLam,dLam,dnu,eNu) !Note: g_tt= -eNu
        use type_mod,only:G,pi,c
        real(8),intent(in)::r,p,rho,m,nu
        real(8),intent(out)::eLam,dLam,dnu,eNu
        eLam = (1.d0 - 2.d0*G*m/r/c**2)**(-1.d0)
        dLam = 2.d0*eLam/r*(G/c**2)*(4.d0*pi*r**2*rho - m/r)
        dnu = 2.d0*eLam/r*(G/c**2)*(4.d0*pi*r**2*P/c**2 + m/r)
        eNu = exp(2.d0*nu)
    end subroutine

    subroutine get_tov_metric2(r,p,rho,m,nu,eLam,dLam,dnu,eNu) !Note: g_tt= -eNu^2
        use type_mod,only:G,pi,c
        real(8),intent(in)::r,p,rho,m,nu
        real(8),intent(out)::eLam,dLam,dnu,eNu
        eLam = (1.d0 - 2.d0*G*m/r/c**2)**(-0.5d0)
        dLam = eLam**2/r*(G/c**2)*(4.d0*pi*r**2*rho - m/r)
        dnu = eLam**2/r*(G/c**2)*(4.d0*pi*r**2*P/c**2 + m/r)
        eNu = exp(nu)
    end subroutine

    subroutine get_bg_r(r,p,rho,m,nu,ga,A,isave)
        use nr_mod,only:polint
        real(8),intent(in)::r
        real(8),intent(out)::p,rho,m,nu,ga,A
        integer,intent(inout)::isave
        real(8),allocatable,dimension(:)::rlist,plist,rholist,mlist,nulist,galist,alist
        real(8)::err
        integer::i,j,itol=5

        rlist=bgs%r
        i=locate_pos(rlist,r,itol,isave)

        isave=max(i-itol,1) !save i for next interpolation
        if (i==1) i=i+1
        if (i+2>size(bgs)) i=size(bgs)-2
        do j=1,size(pt_i)
            if (i==pt_i(j)+1 .and. r>=bgs(pt_i(j))%r) i=i+1
            if (i+2>pt_i(j) .and. r<=bgs(pt_i(j))%r) i=pt_i(j)-2
        end do
        rlist=bgs(i-1:i+2)%r; plist=bgs(i-1:i+2)%p; rholist=bgs(i-1:i+2)%rho
        mlist=bgs(i-1:i+2)%m; nulist=bgs(i-1:i+2)%nu; galist=bgs(i-1:i+2)%ga; alist=bgs(i-1:i+2)%A
        call polint(rlist,plist,r,p,err)
        call polint(rlist,rholist,r,rho,err)
        call polint(rlist,mlist,r,m,err)
        call polint(rlist,nulist,r,nu,err)
        call polint(rlist,galist,r,ga,err)
        call polint(rlist,alist,r,A,err)
    end subroutine

    subroutine get_bg_solid_r(r,p,rho,m,nu,ga,A,mu,chic,isave)
        use nr_mod,only:polint
        real(8),intent(in)::r
        real(8),intent(out)::p,rho,m,nu,ga,A,mu,chic
        integer,intent(inout)::isave
        real(8),allocatable,dimension(:)::rlist,plist,rholist,mlist,nulist,mulist,chiclist,galist,alist
        real(8)::err
        integer::i,j,itol=5

        rlist=bgs%r
        i=locate_pos(rlist,r,itol,isave)

        isave=max(i-itol,1) !save i for next interpolation
        if (i==1) i=i+1
        if (i+2>size(bgs)) i=size(bgs)-2
        do j=1,size(pt_i)
            if (i==pt_i(j)+1 .and. r>=bgs(pt_i(j))%r) i=i+1
            if (i+2>pt_i(j) .and. r<=bgs(pt_i(j))%r) i=pt_i(j)-2
        end do
        rlist=bgs(i-1:i+2)%r; plist=bgs(i-1:i+2)%p; rholist=bgs(i-1:i+2)%rho
        mlist=bgs(i-1:i+2)%m; nulist=bgs(i-1:i+2)%nu; mulist=bgs(i-1:i+2)%mu; chiclist=bgs(i-1:i+2)%chic
        galist=bgs(i-1:i+2)%ga; alist=bgs(i-1:i+2)%A
        call polint(rlist,plist,r,p,err)
        call polint(rlist,rholist,r,rho,err)
        call polint(rlist,mlist,r,m,err)
        call polint(rlist,nulist,r,nu,err)
        call polint(rlist,galist,r,ga,err)
        call polint(rlist,alist,r,A,err)
        call polint(rlist,mulist,r,mu,err)
        call polint(rlist,chiclist,r,chic,err)
    end subroutine

    subroutine get_bg_x(t,r,p,rho,m,nu,ga,A,isave)
        use nr_mod,only:polint
        real(8),intent(in)::t
        real(8),intent(out)::r,p,rho,m,nu,ga,A
        integer,intent(inout)::isave
        real(8),allocatable,dimension(:)::xlist,rlist,plist,rholist,mlist,nulist,galist,alist
        real(8)::err
        integer::i,j,itol=5

        if (.not. allocated(gx)) then
            print*,'err: get_bg_x; gx not allocated'
            stop
        endif
        do while (t<gx(isave))
            isave=max(isave-itol,1) !make sure t is within the range of gx
        end do

        do i=isave,size(bgs)-1
            if (t>=gx(i).and.t<=gx(i+1)) exit
        end do
        if (i==size(bgs)) then
            print*,'err:get_bg_x'
            stop
        endif

        isave=max(i-itol,1) !save i for next interpolation
        if (i==1) i=i+1
        if (i+2>size(bgs)) i=size(bgs)-2
        do j=1,size(pt_i)
            if (i==pt_i(j)+1 .and. t>=gx(pt_i(j))) i=i+1
            if (i+2>pt_i(j) .and. t<=gx(pt_i(j))) i=pt_i(j)-2
        end do
        xlist=gx(i-1:i+2); rlist=bgs(i-1:i+2)%r; plist=bgs(i-1:i+2)%p; rholist=bgs(i-1:i+2)%rho
        mlist=bgs(i-1:i+2)%m; nulist=bgs(i-1:i+2)%nu; galist=bgs(i-1:i+2)%ga; alist=bgs(i-1:i+2)%A
        call polint(xlist,rlist,t,r,err)
        call polint(xlist,plist,t,p,err)
        call polint(xlist,rholist,t,rho,err)
        call polint(xlist,mlist,t,m,err)
        call polint(xlist,nulist,t,nu,err)
        call polint(xlist,galist,t,ga,err)
        call polint(xlist,alist,r,A,err)
    end subroutine

    subroutine get_bg_newt_r(r,p,rho,m,ga,A,isave)
        use nr_mod,only:polint
        real(8),intent(in)::r
        real(8),intent(out)::p,rho,m,ga,A
        integer,intent(inout)::isave
        real(8),allocatable,dimension(:)::rlist,plist,rholist,mlist,galist,alist
        real(8)::err
        integer::i,j,itol=5
!        do while (r<bgs(isave)%r)
!            isave=max(isave-itol,1) !make sure t is within the range of bgs%r
!        end do
!
!        do i=isave,size(bgs)-1
!            if (r>=bgs(i)%r.and.r<=bgs(i+1)%r) exit
!        end do
!        if ((bgs(1)%r-r)/bgs(1)%r>torr_ri) then
!            write(*,*)'err:get_bg: r<<r_i'
!            stop
!        end if
!        if ((r-bgs(size(bgs))%r)/bgs(size(bgs))%r>torr_rf) then
!            write(*,*)'err:get_bg: r>>r_f'
!            stop
!        end if
!        if (i==size(bgs)) then
!            write(*,*)'err:get_bg'
!            stop
!        endif
!        if (i==0) then
!            write(*,*)'err:get_bg'
!            stop
!        endif
        rlist=bgs%r
        i=locate_pos(rlist,r,itol,isave)

        isave=max(i-itol,1) !save i for next interpolation
        if (i==1) i=i+1
        if (i+2>size(bgs)) i=size(bgs)-2
        do j=1,size(pt_i)
            if (i==pt_i(j)+1 .and. r>=bgs(pt_i(j))%r) i=i+1
            if (i+2>pt_i(j) .and. r<=bgs(pt_i(j))%r) i=pt_i(j)-2
        end do
        rlist=bgs(i-1:i+2)%r; plist=bgs(i-1:i+2)%p; rholist=bgs(i-1:i+2)%rho
        mlist=bgs(i-1:i+2)%m; galist=bgs(i-1:i+2)%ga; alist=bgs(i-1:i+2)%A
        call polint(rlist,plist,r,p,err)
        call polint(rlist,rholist,r,rho,err)
        call polint(rlist,mlist,r,m,err)
        call polint(rlist,galist,r,ga,err)
        call polint(rlist,alist,r,A,err)
    end subroutine

    subroutine get_x
        if (allocated(gx)) deallocate(gx)
        allocate(gx(size(bgs)))
        gx=log(bgs%r/bgs%p)
    end subroutine

    function locate_pos(t,tin,itol,isave)
        use io_mod,only:assert
        real(8),intent(in),dimension(:)::t
        real(8),intent(in)::tin
        integer,intent(inout)::isave
        integer,intent(in)::itol
        integer::locate_pos
        integer::i

        do while (tin<t(isave))
            isave=max(isave-itol,1) !make sure tin is within the range of t
        end do

        do i=isave,size(t)-1
            if (tin>=t(i).and.tin<=t(i+1)) exit
        end do

!        if ((t(1)-tin)/t(1)>torr_ri) then
!            write(*,*)'err:locate_pos: tin<<r_i'
!            stop
!        end if
!        if ((tin-t(size(t)))/t(size(t))>torr_rf) then
!            write(*,*)'err:locate_pos: tin>>r_f'
!            stop
!        end if
!        if (i==size(t)) then
!            write(*,*)'err:locate_pos'
!            stop
!        endif
!        if (i==0) then
!            write(*,*)'err:locate_pos'
!            stop
!        endif
        call assert((t(1)-tin)/t(1)<=torr_ri,'err: locate_pos: tin<<r_i')
        call assert((tin-t(size(t)))/t(size(t))<=torr_rf,'err: locate_pos: tin>>rf')
!        call assert(i/=0 .and. i/=size(t),'err:locate_pos')
        call assert(i/=0,'err:locate_pos: i=0')

        locate_pos=i

    end function locate_pos
end module
