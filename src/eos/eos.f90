module eos_mod
    !Using EOS table
    !Default format of the EOS table:
    ! - Col 1 index
    ! - Col 2 baryon number density (not used, mostly)
    ! - Col 3 total energy density in cgs unit
    ! - Col 4 pressure in cgs unit
    implicit none
    private::lnint
    public::etab
    public::eos_read,eos_pt,eos_p,eos_rho,eos_ga_r
    public::eos_ga_r_sm
    public::mu_read,eos_mu,eos_chic

    type eos_table
        integer,allocatable,dimension(:)::i
        real(8),allocatable,dimension(:)::p,rho,nb
        real(8),allocatable,dimension(:)::p2,mu,chi
        integer::n
        integer::n2
    end type eos_table

    type(eos_table)::etab

    real(8),allocatable::eostab(:,:) !First entry: Column; Second entry: Row
    real(8),allocatable::mutab(:,:) !First entry: Column; Second entry: Row

    contains

    subroutine eos_read(path)
        use io_mod,only:file_constTable
        character(len=*)::path
        call file_constTable(path,eostab)
        etab%i      = eostab(1,:)
        etab%nb     = eostab(2,:)
        etab%rho    = eostab(3,:)
        etab%p      = eostab(4,:)
        etab%n      = size(eostab,2)
    end subroutine

    function eos_pt(rhoc,errn) !Look for 1st order phase transitions in the EOS table
        real(8),intent(in)::rhoc
        real(8),allocatable::eos_pt(:),temp(:)
        integer,optional,intent(out)::errn
        integer::i,npt
        allocate(temp(etab%n))
        npt=0
        do i=1,etab%n-1
            if (etab%p(i)==etab%p(i+1).and.rhoc>=etab%rho(i+1)) then
                npt=npt+1
                temp(npt)=etab%p(i)
            endif
            if (etab%p(i)==etab%p(i+1).and.rhoc>etab%rho(i).and.rhoc<etab%rho(i+1)) then
                print*,'err: eos_mod eos_pt: rhoc lies within density gap'
                if (present(errn)) then
                    allocate(eos_pt(0))
                    errn=1
                    return
                else
                    stop
                end if
            endif
        end do
        allocate(eos_pt(npt))
        if (npt/=0) eos_pt(1:npt)=temp(npt:1:-1) !reverse order
        if (present(errn)) errn=0
    end function

    function eos_p(rho)
        real(8),intent(in)::rho
        real(8)::eos_p
        integer::i
        do i=1,etab%n-1
            if (etab%rho(i)<rho .and. etab%rho(i+1)>=rho) exit
        enddo
        if (i==etab%n) then
            print*,'err:eos_p'
            stop
        endif
        eos_p= lnint([etab%rho(i),etab%rho(i+1)],[etab%p(i),etab%p(i+1)],rho)
    end function eos_p

    function eos_rho(p)
        real(8),intent(in)::p
        real(8)::eos_rho
        integer::i
        do i=1,etab%n-1
            if (etab%p(i)<p .and. etab%p(i+1)>=p) exit
        enddo
        if (i==etab%n) then
            print*,'err:eos_rho'
            stop
        endif
        eos_rho= lnint([etab%p(i),etab%p(i+1)],[etab%rho(i),etab%rho(i+1)],p)
    end function eos_rho

    function eos_ga_r(p)
    use type_mod,only:c
        real(8),intent(in)::p
        real(8)::eos_ga_r
        integer::i
        do i=1,etab%n-1
            if (etab%p(i)<p .and. etab%p(i+1)>=p) exit
        enddo
        if (i==etab%n) then
            print*,'err:eos_ga_r'
            stop
        endif
        eos_ga_r= (log(etab%p(i+1))-log(etab%p(i))) &
        /(log(etab%rho(i+1))-log(etab%rho(i)))*(1.d0+p/eos_rho(p)/c**2)
    end function eos_ga_r

    function eos_ga_r_sm(p) !use logistic function to smoothen gamma of the EOS table
    use type_mod,only:c
        real(8),intent(in)::p
        real(8)::eos_ga_r_sm
        real(8),allocatable,dimension(:)::pm
        real(8)::p1,p2,p3
        real(8)::rho1,rho2,rho3
        real(8)::ga1,ga2,ga_p
        real(8)::width
        integer::i
        allocate(pm(etab%n-1))
        do i=1,etab%n-1
            pm(i)=(etab%p(i)+etab%p(i+1))/2
        end do
        do i=1,etab%n-2
            if (pm(i)<p .and. pm(i+1)>=p) exit
        enddo
        if (i==etab%n) then
            print*,'err: eos_ga_r_sm'
            stop
        endif
        p1=etab%p(i);p2=etab%p(i+1);p3=etab%p(i+2)
        rho1=etab%rho(i);rho2=etab%rho(i+1);rho3=etab%rho(i+2)
        ga1=(log(p2)-log(p1))/(log(rho2)-log(rho1))
        ga2=(log(p3)-log(p2))/(log(rho3)-log(rho2))
        if (p2==p3) then
            ga_p=ga1
        elseif (p1==p2) then
            ga_p=ga2
        else
!            width=0.01d0*(p3-p1)
            width=0.005d0*(p3-p1)
            ga_p=ga1+ (ga2-ga1)/(1.d0+exp(-(p-p2)/width)) !logistic function: 1/(1+e^{-(x-x0)/width})
        end if
        eos_ga_r_sm=ga_p*(1.d0+p/eos_rho(p)/c**2)
!eos_ga_r_sm=ga_p*(1.d0+p/eos_rho(p)/c**2)*1.1
    end function eos_ga_r_sm

    !Interpolate shear modulus

    subroutine mu_read(path)
    use io_mod,only:file_constTable
        character(len=*)::path
        call file_constTable(path,mutab)
        etab%p2=mutab(1,:)
        etab%mu=mutab(2,:)
        etab%chi=mutab(3,:)
        etab%n2=size(mutab,2)
    end subroutine

    impure elemental function eos_mu(p)
        real(8),intent(in)::p
        real(8)::eos_mu
        integer::i
        do i=1,etab%n2-1
            if (etab%p2(i)<p .and. etab%p2(i+1)>=p) exit
        enddo
        if (i==etab%n2) then
            print*,'err:eos_mu'
            stop
        endif
        if (etab%mu(i)==0 .and.etab%mu(i+1)==0) then
            eos_mu=0
        elseif (etab%mu(i)==0) then
            eos_mu=etab%mu(i+1)
        else
            eos_mu=lnint([etab%p2(i),etab%p2(i+1)],[etab%mu(i),etab%mu(i+1)],p)
        end if
    end function eos_mu

    impure elemental function eos_chic(p)
        real(8),intent(in)::p
        real(8)::eos_chic
        integer::i
        do i=1,etab%n2-1
            if (etab%p2(i)<p .and. etab%p2(i+1)>=p) exit
        enddo
        if (i==etab%n2) then
            print*,'err:eos_chic'
            stop
        endif
        if (etab%chi(i)==0 .and.etab%chi(i+1)==0) then
            eos_chic=0
        elseif (etab%chi(i)==0) then
            eos_chic=etab%chi(i+1)
        else
            eos_chic=lnint([etab%p2(i),etab%p2(i+1)],[etab%chi(i),etab%chi(i+1)],p)
        end if
    end function eos_chic

    function lnint(xa,ya,x) !linear log interpolation
        real(8),intent(in)::xa(2),ya(2),x
        real(8)::lnint
        lnint=exp((log(ya(2))-log(ya(1)))/(log(xa(2))-log(xa(1)))*(log(x)-log(xa(1)))+log(ya(1)))
    end function lnint
end module
