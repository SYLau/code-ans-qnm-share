module odeSolver_mod
    implicit none
    private
    public::ep,hp,tp,xp
    public::rk45ad
    public::rk4g
    public::rk4ug

    !precision
    integer, parameter::dp=kind(1.d0)
    integer, parameter::qp=kind(1.q0)
    integer, parameter::wp=dp

    real(wp), dimension(:), allocatable :: ep
    real(wp), dimension(:), allocatable :: hp
    real(wp), dimension(:), allocatable :: tp
    real(wp), dimension(:,:), allocatable :: xp

    interface reallocate_a
        module procedure reallocate_rv_a, reallocate_rm_a
    endinterface reallocate_a

    interface
        subroutine input_deriv(t,x,f)
            import wp
            real(wp), intent(in):: t, x(:)
            real(wp), intent(out):: f(:)
        endsubroutine input_deriv

        function term_cond(t,x)
            import wp
            real(wp), intent(in):: t, x(:)
            real(wp):: term_cond
        endfunction term_cond
    endinterface

contains

!===================================================================
! Example block
!===================================================================
!    solve_ode: block
!        use odeSolver_mod,only:rk45ad
!        real(wp)::t,tb
!        real(wp),dimension(:),allocatable::x
!        real(wp)::h
!        real(wp)::emax
!        real(wp)::hmin
!        integer::itmax
!        integer::iflag
!
!        h=!set_hi
!        emax=!set_emax
!        hmin=!set_hmin
!        itmax=!set_itmax
!
!        t=!ti
!        tb=!tf
!        x=!xi
!        call rk45ad(fcn,t,x,h,tb,emax,saveData=.true.,hmin=hmin)
!    end block solve_ode


    !=================================================================================================
    ! RK45 Adaptive using Cash Karp method. Modified from numerical recipes
    !=================================================================================================
    subroutine rk45ad(fcn,t,x,h,tb,emax,saveData,atol,itmax,iflag,xScaleIn,hmin,dtSaveIn,iwarn &
        ,tcond,openbd,tsign)
        ! Modify xScale for different ODEs
        ! n: no of equations; fcn: subroutine containing functions of dx/dt; t: independent variable; x: dependent variables
        ! h: initial guess of step size; itmax: maximum iteration step to increase step size
        ! emin: minimum relative error; emax: maximum relative error; iflag: flag 0 :-> integration finished, flag 1 :-> too many iterations
        ! Optionals: saveData: logical, use RAM to save data or not
        ! Optionals: xScaleIn: fixed scale for relative error, hmin: h automatically set to hmin if abs(h)<abs(hmin) [warning msg is given]
        ! Optionals: iwarn: default = 1. otherwise, suppress warnings
        ! Optionals: tcond: determine terminate conditions. Specified by tcond(t, x(t)) = 0. The zero is determined with a Secant method (1e-8 acc).
        !            openbd: open boundary problem. Whether to disable using tb to terminate.
        implicit none
        procedure(input_deriv)::fcn
        real(wp),intent(in):: tb,emax
        real(wp), intent(inout):: t,x(:),h
        logical, optional, intent(in):: saveData
        real(wp),dimension(size(x)), optional, intent(in):: atol
        integer, optional, intent(in):: itmax
        integer, optional, intent(out):: iflag
        real(wp), optional, intent(in):: xScaleIn(:),hmin
        real(wp), optional, intent(in):: dtSaveIn
        integer,optional,intent(in)::iwarn
        procedure(term_cond),optional:: tcond
        logical,optional::openbd
        integer,optional,intent(in)::tsign
!        real(wp) ::delta=0.5d-14
        logical::savedt
        real(wp),dimension(size(x))::at !absolute tolerance
        real(wp),dimension(size(x))::tol !tolerance
        integer:: itlimit
        integer:: ind !flag
        real(wp),dimension(size(x))::yerr, f
        real(wp):: tsave, dtSave
        real(wp),dimension(size(x))::xsave
        real(wp),dimension(size(x))::xScale
        real(wp)::d
        real(wp)::e
        real(wp)::tcond_old, t_old, t_end, t_acc
        real(wp),dimension(size(x))::x_old
        integer :: n, k
        real(wp), parameter:: safety=0.9_wp, pGrow= -0.2_wp, pShrink= -0.25_wp, errcon=1.89e-4_wp, tiny= 1.e-30_wp !The value errcon equals (5/safety)**(1/pGrow), requires scale-up factor to be at most 5
        logical,save:: warn=.true.
        integer::a_iwarn
        integer::term_sign
        logical::opbd
        integer::dataSize
        real(wp):: tsas

        savedt = .true.
        if (present(saveData)) savedt = saveData

        if (present(xScaleIn)) xScale = xScaleIn

        at = 0
        if (present(atol)) at = atol

        if (present(tcond)) t_old = t !save the previous t, used only when determining tcond = 0
        if (present(tcond)) tcond_old = tcond(t,x)*1.1 !save the previous tcond, used only when determining tcond = 0
        if (present(tcond)) t_acc = 1.e-8_wp !the root-finding method accuracy when solving for tcond = 0

        opbd = .false.
        if (present(openbd)) opbd = openbd ! open boundary

        itlimit = 100
        if (present(itmax)) itlimit = itmax

        a_iwarn = 1
        if (present(iwarn)) a_iwarn = iwarn

        term_sign = 0
        if (present(tsign)) term_sign = tsign

        n = size(x)
        ind = 1 ! ind 0  :-> integration finished
        k = 0
        dataSize = 0
        dtSave = 0._wp
        if (present(dtSaveIn)) dtSave = dtSaveIn
        tsas = t-2._wp*dtSave
        if (allocated(hp)) deallocate(ep)
        if (allocated(hp)) deallocate(hp)
        if (allocated(tp)) deallocate(tp)
        if (allocated(xp)) deallocate(xp)
        allocate(tp(256))
        allocate(hp,ep,source=tp)
        allocate(xp(size(x),size(tp)))

        ! Save first step
        yerr = 0
        if (.not.present(xScaleIn)) then
            call fcn(t, x, f)
            xScale = abs(x)+abs(h*f) + tiny !modify this line to give a custom xScale
        endif
        if (savedt) call save_a_step

        do

            k = k + 1
            if (k > itlimit) exit ! Terminate with error

            if (present(tcond)) then
                if (tcond(t,x)*tcond_old<0) then
!                    t_end = rtsec(terminate_cond,t_old,t,t_acc)
                    t_end = rtbis(terminate_cond,t_old,t,t_acc,term_sign)

                    if (abs(t_end -t_old)<= abs(tb -t_old) .or. opbd) then
                        ind = 0 !signals integration finished
                        t = t_old
                        x = x_old
                        call rkck(fcn,t,x,t_end-t_old,yerr) !integrate a single step using Runge-Kutta-Cash-Karp method

                        if (savedt) call save_a_step
                        if (savedt) call truncate_data
                        exit
                    end if
                end if
            end if

            if (.not. (opbd.and.present(tcond))) then
                d = abs(tb - t)
                if(d <= abs(h))  then !integration range smaller than stepsize  :-> integrate for one more step
                    ind = 0
!                if(d <= delta*max(abs(tb),abs(t))) then
!                    !if (savedt) call save_a_step !no increase in step, simply terminate
!                    if (savedt) call truncate_data
!                    exit    !integration range smaller than bounds * delta :-> integration finished
!                endif
                    h = sign(1._wp,h)*d     !make the last stepsize = d
                end if
            end if

            xsave = x
            tsave = t
            t_old = t
            x_old = x
            if (present(tcond)) tcond_old = tcond(t,x)

            if (.not.present(xScaleIn)) then
                call fcn(t, x, f)
                xScale = abs(x)+abs(h*f) + tiny !modify this line to give a custom xScale
            endif

            call rkck(fcn,t,x,h,yerr) !integrate a single step using Runge-Kutta-Cash-Karp method


            tol = emax*xScale + at
            e = maxval(abs(yerr/tol))

            if (e > 1) ind = 1 !check if e exceeds tol
            if (ind == 0) then !exit and save last step
                if (savedt) call save_a_step
                if (savedt) call truncate_data
                exit
            endif

            if (present(hmin)) then !new
                if (abs(h)<abs(hmin)) then
                    if (warn) then
                        if (a_iwarn == 1) print*, 'warning:rk45ad stepsize smaller than hmin; replacing h by hmin'
                    endif
                    warn = .false.
                    k = 0
                    if (abs(t-tsas) > abs(dtSave) .and. savedt) call save_a_step
                    h=sign(abs(hmin),h)
                    cycle
                endif
            endif

            if (e > 1) then !integrate again with shrinked step size
                h= sign(max(abs(safety*h*(e)**pShrink), abs(h)/10._wp), h)
                x = xsave
                t = tsave
                if(t == t+h) print*,'error occurs at:', t,h
                if(t == t+h) print*,e, yerr, tol
                if(t == t+h) stop 'stepsize underflow rk45ad. Probably stiff ODE. Or bad local behaviours in ODE'

            else

                k = 0
                if (abs(t-tsas) > abs(dtSave) .and. savedt) call save_a_step
                if (e > errcon) then !check if the scale-up process is too fast
                    h= safety*h*(e)**pGrow !next step size is scaled up
                else
                    h= 5._wp*h !next step size is scaled up by at most 5
                endif

            endif
            !if (present(hmin)) then
            !  if (abs(h)<abs(hmin)) stop 'err:rk45ad stepsize smaller than hmin'
            !endif

        enddo

        warn = .true.
        if (present(iflag)) iflag = ind
        if (ind /= 0) stop 'ode_solver: err: exceeds maximum iterations'

    contains

        subroutine save_a_step
            dataSize = dataSize +1
            if (size(tp)<dataSize) then
                tp=reallocate_a(tp,2*size(tp))
                hp=reallocate_a(hp,2*size(hp))
                ep=reallocate_a(ep,2*size(ep))
                xp=reallocate_a(xp,size(xp,1),2*size(tp))
            endif
            tp(dataSize)=t
            hp(dataSize)=h
            ep(dataSize)=maxval(abs(yerr/xScale))
            xp(:,dataSize)=x(:)
            tsas = t
        end subroutine

        subroutine truncate_data
            tp=reallocate_a(tp,dataSize)
            hp=reallocate_a(hp,dataSize)
            ep=reallocate_a(ep,dataSize)
            xp=reallocate_a(xp,size(xp,1),dataSize)
        end subroutine

        function terminate_cond(tin) result(y)
            real(wp), intent(in) :: tin
            real(wp) :: y
            real(wp)::step
            real(wp)::tstart
            real(wp), dimension(size(x)) ::xt, yerr

            tstart = t_old
            step = tin - tstart
            xt = x_old
            call rkck(fcn,tstart,xt,step,yerr) !integrate a single step using Runge-Kutta-Cash-Karp method

            y = tcond(tin,xt)
        end function terminate_cond

    end subroutine rk45ad

    ! Modified from numerical recipes
    subroutine rkck(derivs,x,y,h,yerr)
        implicit none
        real(wp), dimension(:), intent(inout) :: y
        real(wp), intent(inout) :: x
        real(wp), intent(in) :: h
        real(wp), dimension(:), intent(out) :: yerr
        interface
            subroutine derivs(x,y,dydx)
            import wp
            implicit none
            real(wp), intent(in) :: x
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface
        integer :: ndum
        real(wp), dimension(size(y)) :: dydx
        real(wp), dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
        real(wp), parameter :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
            A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,&
            B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,&
            B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,&
            B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,&
            B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,&
            B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,&
            C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,&
            C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,&
            DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,&
            DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp

        call derivs(x,y,dydx)

        ytemp=y+B21*h*dydx
        call derivs(x+A2*h,ytemp,ak2)
        ytemp=y+h*(B31*dydx+B32*ak2)
        call derivs(x+A3*h,ytemp,ak3)
        ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs(x+A4*h,ytemp,ak4)
        ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs(x+A5*h,ytemp,ak5)
        ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs(x+A6*h,ytemp,ak6)
        x=x+h
        y=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
        yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
    end subroutine rkck

    ! Replaced by Numerical Recipes rkck
    subroutine rkck45(fcn,t,x,h,yerr)   !Runge-Kutta-Cash-Karp method
        implicit none
        real(wp), dimension(:), intent(inout) :: x
        real(wp), intent(inout) :: t
        real(wp), intent(in) :: h
        real(wp), dimension(:), intent(out) :: yerr
        integer :: n
        real(wp), dimension(1:size(x)):: xi, f, f1,f2,f3,f4,f5,f6, x5
        real(wp):: c20,c21,c30,c31,c32,c40,c41,c42,c43,c51,c52,c53,c54,c60,c61,c62,c63,c64,c65,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6
        interface deriv
            subroutine fcn(t,x,f)
                import wp
                real(wp), intent(in):: t, x(:)
                real(wp), intent(out):: f(:)
            endsubroutine fcn
        endinterface deriv

        c20=0.2_wp; c21=0.2_wp; c30=0.3_wp; c31=0.075_wp; c32=0.225_wp
        c40=3._wp/5._wp; c41=3._wp/10._wp; c42=-9._wp/10._wp; c43=6._wp/5._wp
        c51=-11._wp/54._wp; c52=5._wp/2._wp; c53=-70._wp/27._wp; c54=35._wp/27._wp
        c60=0.875_wp; c61=1631._wp/55296; c62=175._wp/512; c63=575._wp/13824; c64=44275._wp/110592; c65=253._wp/4096
        a1=2825._wp/27648._wp; a2=0._wp; a3=18575._wp/48384._wp; a4=13525._wp/55296._wp; a5=277._wp/14336._wp; a6=0.25_wp
        b1=37._wp/378._wp; b2=0._wp; b3=250._wp/621._wp; b4=125._wp/594._wp; b5=0._wp; b6=512._wp/1771._wp
        n = size(x)
        xi(1:n) = x(1:n)
        call fcn(t, xi, f)
        f1(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c21*f1(1:n)
        call fcn(t+ c20*h, xi, f)
        f2(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c31*f1(1:n) + c32*f2(1:n)
        call fcn(t+ c30*h, xi, f)
        f3(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c41*f1(1:n) + c42*f2(1:n) + c43*f3(1:n)
        call fcn(t+ c40*h, xi, f)
        f4(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c51*f1(1:n) + c52*f2(1:n) + c53*f3(1:n) + c54*f4(1:n)
        call fcn(t+h, xi, f)
        f5(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c61*f1(1:n) + c62*f2(1:n) + c63*f3(1:n) + c64*f4(1:n) + c65*f5(1:n)
        call fcn(t+ c60*h, xi, f)
        f6(1:n) = h*f(1:n)

        x5(1:n) = x(1:n) + b1*f1(1:n) + b3*f3(1:n) + b4*f4(1:n) + b5*f5(1:n) + b6*f6(1:n)
        x(1:n) = x(1:n) + a1*f1(1:n) + a3*f3(1:n) + a4*f4(1:n) + a5*f5(1:n)  +  a6*f6(1:n)
        t = t + h
        yerr(1:n) = x(1:n) - x5(1:n)
    end subroutine rkck45


    !=================================================================================================
    ! RK4 global error control
    !=================================================================================================
    subroutine rk4g(fcn,t,x,tb,saveData,nlist,nmin,itmax,rtol)
        real(wp), dimension(:),intent(inout)::x
        real(wp),intent(inout)::t
        real(wp), intent(in) ::tb
        logical,intent(in),optional::saveData
        integer,intent(in),optional,dimension(:)::nlist
        integer,intent(in),optional::nmin
        integer,intent(in),optional::itmax
        real(wp),optional,intent(in)::rtol
        integer,allocatable,dimension(:)::steps
        integer::j
        real(wp)::t0
        real(wp), dimension(size(x))::x0
        real(wp), dimension(size(x))::xold
        real(wp), dimension(size(x))::tol

        logical::savedt
        real(wp)::eps
        integer::itlimit
        integer::nstep
        real(wp),parameter::tiny = 1.e-30_wp
        interface deriv
            subroutine fcn(t,x,f)
                import wp
                real(wp), intent(in):: t, x(:)
                real(wp), intent(out):: f(:)
            endsubroutine fcn
        endinterface deriv

        savedt=.true.
        if (present(saveData)) savedt = saveData
        nstep=500
        if (present(nmin)) nstep = nmin
        itlimit=15
        if (present(itmax)) itlimit = itmax
        eps=1.e-4_wp
        if (present(rtol)) eps = rtol
        steps=[(j*nstep,j=1,itlimit)]
        if (present(nlist)) then
            steps=nlist
            if (size(steps)<2) stop 'err: rk4g nlist < 2'
        end if

        t0=t
        x0=x

        xold=0
        do j=1,size(steps)
            t=t0
            x=x0
            call rk4ug(fcn,t,x,steps(j),tb,savedt)

            tol = eps*max(abs(x),abs(xold)) + tiny
            if (maxval(abs((x-xold))/tol) <= 1._wp ) return
            xold=x
        end do

        write(*,*) 'rk4g: too many steps'
        stop
    end subroutine rk4g

    !=================================================================================================
    ! RK4 uniform grid
    !=================================================================================================
    subroutine rk4ug(fcn,t,x,steps,tb,saveData)
        implicit none
        real(wp), dimension(:),intent(inout)::x
        real(wp),intent(inout)::t
        real(wp), intent(in) ::tb
        integer,intent(in)::steps
        logical,intent(in),optional::saveData
        real(wp)::ti,h
        integer::i
        logical::savedt
        interface deriv
            subroutine fcn(t,x,f)
                import wp
                real(wp), intent(in):: t, x(:)
                real(wp), intent(out):: f(:)
            endsubroutine fcn
        endinterface deriv

        savedt=.true.
        if (present(saveData)) savedt = saveData

        if (allocated(tp)) deallocate(tp)
        if (allocated(xp)) deallocate(xp)
        allocate(tp(256))
        allocate(xp(size(x),size(tp)))
        ti=t
        h=(tb-ti)/steps

        do i=1,steps
            t=ti+(i-1)*h
            call rk4(fcn,t,x,h)
            if (savedt) call save_a_step
        end do
        if (savedt) call truncate_data

        contains
        subroutine save_a_step
            if (size(tp)<i) then
                tp=reallocate_a(tp,2*size(tp))
                xp=reallocate_a(xp,size(xp,1),2*size(tp))
            endif
            tp(i)=t
            xp(:,i)=x(:)
        end subroutine
        subroutine truncate_data
            tp=reallocate_a(tp,steps)
            xp=reallocate_a(xp,size(xp,1),steps)
        end subroutine
    end subroutine

    subroutine rk4(fcn,t,x,h) !Return one step of RK4
        implicit none
        real(wp), dimension(:), intent(inout) :: x
        real(wp), intent(in) :: t,h
        interface deriv
            subroutine fcn(t,x,f)
                import wp
                real(wp), intent(in):: t, x(:)
                real(wp), intent(out):: f(:)
            endsubroutine fcn
        endinterface deriv
        real(wp) :: h6,hh,th
        real(wp), dimension(size(x)) :: dxdt,dxm,dxt,xt
        hh=h*0.5_wp
        h6=h/6._wp
        th=t+hh
        call fcn(t,x,dxdt)
        xt=x+hh*dxdt
        call fcn(th,xt,dxt)
        xt=x+hh*dxt
        call fcn(th,xt,dxm)
        xt=x+h*dxm
        dxm=dxt+dxm
        call fcn(t+h,xt,dxt)
        x=x+h6*(dxdt+dxt+2._wp*dxm)
    endsubroutine rk4

    !=================================================================================================
    ! Utilities
    !=================================================================================================
    function reallocate_rv_a(p,n)
        real(wp), dimension(:), allocatable :: p, reallocate_rv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_rv_a(n),stat=ierr)
        if (ierr /= 0) stop 'odeSolver_mod: reallocate_rv_a: problem in attempt to allocate memory'
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_rv_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_rv_a
    function reallocate_rm_a(p,n,m)
        real(wp), dimension(:,:), allocatable :: p, reallocate_rm_a
        integer, intent(in) :: n,m
        integer :: nold,mold,ierr
        allocate(reallocate_rm_a(n,m),stat=ierr)
        if (ierr /= 0) stop 'odeSolver_mod: reallocate_rm_a: problem in attempt to allocate memory'
        if (.not. allocated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm_a(1:min(nold,n),1:min(mold,m))=p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    endfunction reallocate_rm_a
    !------------------------------------------------------------------------------
    ! Local Utility
    !------------------------------------------------------------------------------
    ! Modified from numerical recipes
    function rtbis(func,x1,x2,xacc,term_sign)
        implicit none
        real(wp), intent(in) :: x1,x2,xacc
        integer,intent(in)::term_sign
        real(wp) :: rtbis
        interface
            function func(x)
                import wp
                implicit none
                real(wp), intent(in) :: x
                real(wp) :: func
            end function func
        end interface
        integer, parameter :: MAXIT=150
        integer :: j
        real(wp) :: dx,f,fmid,xmid
        fmid=func(x2)
        f=func(x1)

        if (f*fmid >= 0.0) stop 'rtbis: root must be bracketed'
        if (f < 0.0) then
            rtbis=x1
            dx=x2-x1
        else
            rtbis=x2
            dx=x1-x2
        end if
        do j=1,MAXIT
            dx=dx*0.5_dp
            xmid=rtbis+dx
            fmid=func(xmid)
            if (fmid <= 0.0) rtbis=xmid
!            if (abs(dx) < xacc .or. fmid == 0.0) return
            if (abs(dx) < xacc .or. fmid == 0.0) then
                if (term_sign == 0) then
                    return
                else
                    if (fmid*sign(1,term_sign)>=0.) then
                        rtbis = xmid
                        return
                    else
                        rtbis = xmid+dx
                        return
                    endif
                endif
            endif
        end do
        stop 'ode_solver: rtbis: too many bisections'
    end function rtbis

    !  Modified from secant method from numerical recipes
    function rtsec(func,x1,x2,xacc)
        real(wp), intent(in) :: x1,x2,xacc
        real(wp) :: rtsec
        integer, parameter :: MAXIT=30
        integer :: j
        real(wp) :: dx,f,fl,xl
        real(wp)::fsave,xsave
        interface
            function func(x)
                import wp
                real(wp), intent(in) :: x
                real(wp) :: func
            end function func
        end interface
        fl=func(x1)
        f=func(x2)
        if (abs(fl) < abs(f)) then
            rtsec=x1
            xl=x2
            call swap(fl,f)
        else
            xl=x1
            rtsec=x2
        end if
        do j=1,MAXIT
            dx=(xl-rtsec)*f/(f-fl)
            xsave = xl
            fsave = fl
            xl=rtsec
            fl=f
            rtsec=rtsec+dx
            f=func(rtsec)
            if (abs(dx) < xacc .or. f == 0.0) return
        end do
        stop 'ode_solver: rtsec: exceed maximum iterations'
    end function rtsec
    subroutine swap(a,b)
        real(wp), intent(inout) :: a,b
        real(wp) :: dum
        dum=a
        a=b
        b=dum
    end subroutine swap

endmodule odeSolver_mod

