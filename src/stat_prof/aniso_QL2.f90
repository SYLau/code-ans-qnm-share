module anisotropy_QL_2_mod
    use type_mod,only:wp
    use anisotropy_type_mod,only: set_par, sigma_func, dsigma_func, sigma_coef
    implicit none
    private
    public:: aniso_QL2
    !QL model ref: https://iopscience.iop.org/article/10.1088/0264-9381/28/2/025009
    !Modified QL model to satisfy the requirement in full GR



!    interface
!        subroutine set_par(x)
!            import wp
!            real(wp),intent(in)::x
!        end subroutine set_par
!
!        function sigma_func(p,rho,m,r) result(res)
!            import wp
!            real(wp),intent(in)::p,rho,m,r
!            real(wp)::res
!        end function sigma_func
!
!        function sigma_coef(p,rho) result(res)
!            import wp
!            real(wp),intent(in)::p,rho
!            real(wp)::res
!        end function sigma_coef
!
!    end interface

    real(wp)::lam = 0

    type aniso_QL2
        procedure(set_par), pointer, nopass :: set              => set_lam

        procedure(sigma_func), pointer, nopass :: s             => sigma_QL2
        procedure(dsigma_func), pointer, nopass :: dsp          => dsigma_p_QL2
        procedure(dsigma_func), pointer, nopass :: dsrho        => dsigma_rho_QL2
        procedure(dsigma_func), pointer, nopass :: dsmu         => dsigma_mu_QL2

        procedure(dsigma_func), pointer, nopass :: ddsp         => ddsigma_p_QL2
        procedure(dsigma_func), pointer, nopass :: ddspmu       => ddsigma_p_mu_QL2
        procedure(dsigma_func), pointer, nopass :: dmu          => dmu_r_QL2
        procedure(dsigma_func), pointer, nopass :: ds           => dsigma_r_QL2

        procedure(sigma_coef), pointer, nopass :: s2            => sigma_2_QL2
    end type aniso_QL2


    contains

    subroutine set_lam(lam_in)
        real(wp),intent(in)::lam_in
        lam=lam_in
    end subroutine set_lam

    function sigma_QL2(p,rho,m,r) result(res) !geometrized unit
        real(wp),intent(in)::p,rho,m,r
        real(wp)::mu
        real(wp)::res
        mu=2*m/r
        res = lam*p*mu**2
    end function sigma_QL2

    function dsigma_p_QL2(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        real(wp)::mu
        mu=2*m/r
        res = lam*mu**2
    end function dsigma_p_QL2

    function dsigma_rho_QL2(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res = 0
    end function dsigma_rho_QL2

    function dsigma_mu_QL2(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::mu
        real(wp)::res
            mu=2*m/r
            res = 2*lam*p*mu
    end function dsigma_mu_QL2

    !================================================================================
    ! 2nd derivatives of sigma; won't needed if choosing dP as dependent variable instead of V
    !================================================================================
    function ddsigma_p_QL2(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res = 0
    end function ddsigma_p_QL2

    function ddsigma_p_mu_QL2(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        real(wp)::mu
        mu=2*m/r
        res = lam*mu
    end function ddsigma_p_mu_QL2

    function dmu_r_QL2(p,rho,m,ga,A,r) result(res)
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res=pi*8*rho*r-2*m/r**2
    end function dmu_r_QL2

    function dsigma_r_QL2(p,rho,m,ga,A,r) result(res)
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::mu
        real(wp)::dp
        real(wp)::res
        mu=2*m/r
        dp=-(rho+p)*(m+4*pi*p*r**3)/r**2/(1-2*m/r) -2*sigma_QL2(p,rho,m,r)/r
        res=(lam*mu**2)*dp+(2*lam*p*mu)*dmu_r_QL2(p,rho,m,ga,A,r)
    end function dsigma_r_QL2

    !================================================================================
    ! Expansion coefficient of sigma near r = 0
    !================================================================================
    function sigma_2_QL2(p,rho,ga,A) result(res)
        ! expansion coefficient at r = 0: sigma = sigma_0 + sigma_2 r^2 + ...
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,ga,A
        real(wp)::res
        res = 0
    end function sigma_2_QL2

end module anisotropy_QL_2_mod
