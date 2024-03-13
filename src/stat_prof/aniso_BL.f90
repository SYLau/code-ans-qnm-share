module anisotropy_BL_mod
    use type_mod,only:wp
    use anisotropy_type_mod,only: set_par, sigma_func, dsigma_func, sigma_coef
    implicit none
    private
    public:: aniso_BL
    !BL model ref: Bowers and Liang 1974ApJ. . .188. .657B
    !Equations (3.1) (3.2) and n = 2.

    real(wp)::lam = 0

    type aniso_BL
        procedure(set_par), pointer, nopass :: set              => set_lam

        procedure(sigma_func), pointer, nopass :: s             => sigma_BL
        procedure(dsigma_func), pointer, nopass :: dsp          => dsigma_p_BL
        procedure(dsigma_func), pointer, nopass :: dsrho        => dsigma_rho_BL
        procedure(dsigma_func), pointer, nopass :: dsmu         => dsigma_mu_BL

!        procedure(dsigma_func), pointer, nopass :: ddsp         => ddsigma_p_BL
!        procedure(dsigma_func), pointer, nopass :: ddspmu       => ddsigma_p_mu_BL
        procedure(dsigma_func), pointer, nopass :: dmu          => dmu_r_BL
        procedure(dsigma_func), pointer, nopass :: ds           => dsigma_r_BL

        procedure(sigma_coef), pointer, nopass :: s2            => sigma_2_BL
    end type aniso_BL


    contains

    subroutine set_lam(lam_in)
        real(wp),intent(in)::lam_in
        lam=lam_in
    end subroutine set_lam

    function sigma_BL(p,rho,m,r) result(res) !geometrized unit
        real(wp),intent(in)::p,rho,m,r
        real(wp)::mu
        real(wp)::res
        mu=2*m/r
        res = lam*(rho+p)*(rho+3*p)/(1-mu)*r**2
    end function sigma_BL

    function dsigma_p_BL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        real(wp)::mu
        mu=2*m/r
        res = 2*lam*(2*rho+3*p)/(1-mu)*r**2
    end function dsigma_p_BL

    function dsigma_rho_BL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        real(wp)::mu
        mu=2*m/r
        res = 2*lam*(rho+2*p)/(1-mu)*r**2
    end function dsigma_rho_BL

    function dsigma_mu_BL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::mu
        real(wp)::res
            mu=2*m/r
            res = -lam*(rho+p)*(rho+3*p)/(1-mu)**2*r**2
    end function dsigma_mu_BL

    !================================================================================
    ! d/dr
    !================================================================================

    function dmu_r_BL(p,rho,m,ga,A,r) result(res)
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res=pi*8*rho*r-2*m/r**2
    end function dmu_r_BL

    function dsigma_r_BL(p,rho,m,ga,A,r) result(res)
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::mu
        real(wp)::cs2
        real(wp)::dp
        real(wp)::drho
        real(wp)::res
        mu=2*m/r
        cs2=ga*p/(rho+p)
        dp=-(rho+p)*(m+4*pi*p*r**3)/r**2/(1-2*m/r) -2*sigma_BL(p,rho,m,r)/r
        drho=A*r +dp/cs2
        res= dsigma_rho_BL(p,rho,m,ga,A,r)*drho +dsigma_p_BL(p,rho,m,ga,A,r)*dp &
        +dsigma_mu_BL(p,rho,m,ga,A,r)*dmu_r_BL(p,rho,m,ga,A,r) + 2*lam*(rho+p)*(rho+3*p)/(1-mu)*r
    end function dsigma_r_BL

    !================================================================================
    ! Expansion coefficient of sigma near r = 0
    !================================================================================
    function sigma_2_BL(p,rho,ga,A) result(res)
        ! expansion coefficient at r = 0: sigma = sigma_0 + sigma_2 r^2 + ...
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,ga,A
        real(wp)::res
        res = lam*(rho+p)*(rho+3*p)
    end function sigma_2_BL

end module anisotropy_BL_mod
