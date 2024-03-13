module odepdr_mod
    use type_mod,only: wp
    use odeSolver_mod, only: tp, xp                                 !Local modifications
    implicit none
    private
    public::odelsoda
!    public::tp, xp                                                 !Local modifications

!    integer, parameter:: dp = kind(1.d0)                           !Local modifications
!    integer, parameter:: wp = dp                                   !Local modifications

!    real(wp), dimension(:), allocatable :: tp                      !Local modifications
!    real(wp), dimension(:,:), allocatable :: xp                    !Local modifications

    interface
        subroutine input_func(t,y,ydot)
            import wp
            real(wp),intent(in):: t
            real(wp), dimension(:),intent(in):: y
            real(wp), dimension(:),intent(out):: ydot
        end subroutine input_func

        subroutine input_jac(t, y, pd)
            import :: wp
            real(wp), intent(in) :: t
            real(wp), intent(inout) :: y(:)
            real(wp), intent(out) :: pd(:,:)
        end subroutine
    end interface

    interface
        subroutine lsoda_f(neq, t, y, ydot)
            import :: wp
            integer, intent(in):: neq
            real(wp), intent(in):: t
            real(wp), dimension(neq), intent(in):: y
            real(wp), dimension(neq), intent(out):: ydot
        end subroutine lsoda_f

        subroutine lsoda_jac(neq, t, y, ml, mu, pd, nrowpd)
            import :: wp
            integer, intent(in):: neq, ml, mu, nrowpd
            real(wp), intent(in):: t
            real(wp), dimension(neq), intent(inout):: y
            real(wp), dimension(nrowpd, neq), intent(out):: pd
        end subroutine lsoda_jac
    end interface

    interface
        subroutine dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask, &
                istate, iopt, rwork, lrw, iwork, liw, jac, jt)
            import :: wp,  lsoda_f, lsoda_jac
            implicit none
            procedure(lsoda_f) :: f
            integer, dimension(*), intent(in) :: neq
            real(wp), dimension(*), intent(inout) :: y
            real(wp), intent(inout) :: t
            real(wp), intent(in) :: tout
            integer, intent(in) :: itol
            real(wp), dimension(*), intent(in) :: rtol
            real(wp), dimension(*), intent(in) :: atol
            integer, intent(in) :: itask
            integer, intent(inout) :: istate
            integer, intent(in) :: iopt
            integer, intent(in) :: lrw
            real(wp), intent(inout) :: rwork(lrw)
            integer, intent(in) :: liw
            integer, intent(inout) :: iwork(liw)
            procedure(lsoda_jac) :: jac
            integer, intent(in) :: jt
        end subroutine

    end interface

contains

    subroutine odelsoda(func,t,y,tf,eps,t_len,min_dt,t_list,abstol,jaco)
        procedure(input_func)::func
        real(wp),intent(inout)::t
        real(wp),dimension(:),intent(inout)::y
        real(wp),intent(in)::tf
        real(wp),intent(in)::eps

        integer,optional, intent(in)::t_len
        real(wp),optional, intent(in)::min_dt
        real(wp),dimension(:),optional, intent(in)::t_list
        real(wp),dimension(:),optional, intent(in)::abstol
        procedure(input_jac), optional::jaco

        integer::neq
        integer::itol
        real(wp)::tout
        real(wp)::rtol
        real(wp),allocatable,dimension(:)::atol
        integer::itask
        integer::istate
        integer::iopt
        integer:: lrw
        real(wp),allocatable,dimension(:):: rwork
        integer:: liw
        integer,allocatable,dimension(:):: iwork
        integer:: jt
        procedure(lsoda_jac), pointer:: jac_choice      => null()

        integer::alen
        real(wp),allocatable,dimension(:)::t_eval

        neq = size(y)
        rtol = eps                                                  !relative tolerance

        assign_atol:block
            integer::i
            real(wp),dimension(size(y))::ydot

            itol = 1                                                !tolerance options (scalars atol or not)
            atol=[0._wp]                                            !default absolute tolerance

            do i=1,size(y)
                if (abs(y(i))<=tiny(1._wp)) then
                    itol = 2
                    atol = y*0
                end if
            end do
            call func(t,y,ydot)
            do i=1,size(y)
                if (abs(y(i))<=tiny(1._wp)) then
                    atol(i) = abs(ydot(i)*(tf-t))*1.e-15
                    if (abs(ydot(i))<=tiny(1._wp)) atol(i)= 1._wp   !for a 0 = 0 ODE
                end if
            end do

            if (present(abstol)) then
                atol=abstol
                if (size(abstol)>1) itol = 2
            endif
        end block assign_atol

        itask = 1                                                   !whether to overshoot tout etc
        istate = 1                                                  !inout: initial call
        iopt = 1                                                    !whether to enable optional inputs

        jt = 2                                                      !jacobian input options: 1, full jacobian input 2, full jacobian with finite differencing
        jac_choice => jac_dum
        if (present(jaco)) jt = 1
        if (present(jaco)) jac_choice => jac

        alen = 2
        if (present(t_len))alen = t_len
        if (.not.present(t_len).and.present(min_dt)) &
            alen= ceiling(abs((tf-t)/min_dt))

        t_eval = linspace(t, tf, alen)
        if (present(t_list)) then
            if ((t_list(1)-t)*(t_list(size(t_list))-t)> 0) &
                error stop 'err: t_list outside t and tf'

            deallocate(t_eval)
            allocate(t_eval(2+size(t_list)))
            t_eval(1)=t
            t_eval(2:size(t_eval)-1)=t_list
            t_eval(size(t_eval))=tf
        end if

        rwork_len: block
            integer::lrn,lrs
            if (jt == 1 .or. jt == 2) then
                lrn = 22 + 16*neq
                lrs = 22 + 9*neq + neq**2
            else
                print*, 'jt = ', jt
                error stop 'err: lsoda jt not supported'
            end if
            lrw = max(lrn,lrs)
        end block rwork_len
        allocate(rwork(lrw))
        rwork=0

        liw = 20 + neq
        allocate(iwork(liw))                                        !inout: band diagonal jac, iwork(1/2) are the lower and upper bandwidth
        iwork=0
        iwork(6) = 10000                                             !optional input: maximum steps

        integrate_w_lsoda: block
            integer::i

            if (allocated(xp)) deallocate(xp)
            allocate(xp(size(y),alen))
            xp(:,1) = y
            do i = 2,alen
                tout = t_eval(i)
                call dlsoda(f,[neq],y,t,tout,itol,[rtol],atol,itask,istate, &
                iopt,rwork,lrw,iwork,liw,jac_choice,jt)

                if (istate < 0) then
                    print*,'istate = ', istate
                    error stop 'err: odelsoda dlsoda istate < 0'
                end if

                xp(:,i) = y
            end do
            tp = t_eval

        end block integrate_w_lsoda

    contains

        subroutine f(neq, t, y, ydot)
            integer, intent(in):: neq
            real(wp), intent(in):: t
            real(wp), dimension(neq), intent(in):: y
            real(wp), dimension(neq), intent(out):: ydot

            call func(t,y,ydot)
        end subroutine f

        subroutine jac(neq, t, y, ml, mu, pd, nrowpd)
            integer, intent(in):: neq, ml, mu, nrowpd
            real(wp), intent(in):: t
            real(wp), dimension(neq), intent(inout):: y
            real(wp), dimension(nrowpd, neq), intent(out):: pd

            call jaco(t, y, pd)
        end subroutine jac

        subroutine jac_dum(neq, t, y, ml, mu, pd, nrowpd)
            integer, intent(in):: neq, ml, mu, nrowpd
            real(wp), intent(in):: t
            real(wp), dimension(neq), intent(inout):: y
            real(wp), dimension(nrowpd, neq), intent(out):: pd
            return
        end subroutine jac_dum

    end subroutine odelsoda



    !==================================================================================
    ! Utilities
    !==================================================================================
    function linspace(from, to, n) result(array)
        real(wp), intent(in) :: from, to
        integer:: n
        real(wp), allocatable, dimension(:) :: array
        real(wp) :: range
        integer :: i

        allocate(array(n))
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if

        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end function linspace

!    subroutine linspace(from, to, array)
!        real(wp), intent(in) :: from, to
!        real(wp), intent(out) :: array(:)
!        real(wp) :: range
!        integer :: n, i
!        n = size(array)
!        range = to - from
!
!        if (n == 0) return
!
!        if (n == 1) then
!            array(1) = from
!            return
!        end if
!
!        do i=1, n
!            array(i) = from + range * (i - 1) / (n - 1)
!        end do
!    end subroutine linspace

end module odepdr_mod
