module anisotropy_QL_mod
    use type_mod,only:wp
    use anisotropy_type_mod,only: set_par, sigma_func, dsigma_func, sigma_coef
    implicit none
    private
    public:: aniso_QL
    !QL model ref: https://iopscience.iop.org/article/10.1088/0264-9381/28/2/025009


!    type aniso_QL
!        real(wp), public::lam = 0
!
!        contains
!
!        procedure, public :: TOV_s      => TOV_sigma_QL
!        procedure, public :: s          => sigma_QL
!        procedure, public :: dsp        => dsigma_p_QL
!        procedure, public :: ddsp       => ddsigma_p_QL
!        procedure, public :: ddspmu     => ddsigma_p_mu_QL
!        procedure, public :: dmu        => dmu_r_QL
!    end type aniso_QL

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

    type aniso_QL
        procedure(set_par), pointer, nopass :: set              => set_lam

        procedure(sigma_func), pointer, nopass :: s             => sigma_QL
        procedure(dsigma_func), pointer, nopass :: dsp          => dsigma_p_QL
        procedure(dsigma_func), pointer, nopass :: dsrho        => dsigma_rho_QL
        procedure(dsigma_func), pointer, nopass :: dsmu         => dsigma_mu_QL

        procedure(dsigma_func), pointer, nopass :: ddsp         => ddsigma_p_QL
        procedure(dsigma_func), pointer, nopass :: ddspmu       => ddsigma_p_mu_QL
        procedure(dsigma_func), pointer, nopass :: dmu          => dmu_r_QL
        procedure(dsigma_func), pointer, nopass :: ds           => dsigma_r_QL

        procedure(sigma_coef), pointer, nopass :: s2            => sigma_2_QL
    end type aniso_QL


    contains

    subroutine set_lam(lam_in)
        real(wp),intent(in)::lam_in
        lam=lam_in
    end subroutine set_lam

    function sigma_QL(p,rho,m,r) result(res) !geometrized unit
        real(wp),intent(in)::p,rho,m,r
        real(wp)::mu
        real(wp)::res
        mu=2*m/r
        res = lam*p*mu
    end function sigma_QL

    function dsigma_p_QL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        real(wp)::mu
        mu=2*m/r
        res = lam*mu
!        temp:block
!!            res = lam*mu**2 !old
!        end block temp

    end function dsigma_p_QL

    function dsigma_rho_QL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res = 0
    end function dsigma_rho_QL

    function dsigma_mu_QL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res = lam*p
!        temp:block
!            real(wp)::mu
!            mu=2*m/r
!!            res = 2*lam*p*mu !old
!            res = 0
!        end block temp
    end function dsigma_mu_QL

    !================================================================================
    ! 2nd derivatives of sigma; won't needed if choosing dP as dependent variable instead of V
    !================================================================================
    function ddsigma_p_QL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res = 0
    end function ddsigma_p_QL

    function ddsigma_p_mu_QL(p,rho,m,ga,A,r) result(res)
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res = lam
    end function ddsigma_p_mu_QL

    function dmu_r_QL(p,rho,m,ga,A,r) result(res)
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::res
        res=pi*8*rho*r-2*m/r**2
    end function dmu_r_QL

    function dsigma_r_QL(p,rho,m,ga,A,r) result(res)
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,m,ga,A,r
        real(wp)::mu
        real(wp)::dp
        real(wp)::res
        mu=2*m/r
        dp=-(rho+p)*(m+4*pi*p*r**3)/r**2/(1-2*m/r) -2*sigma_QL(p,rho,m,r)/r
        res=(lam*mu)*dp+(lam*p)*dmu_r_QL(p,rho,m,ga,A,r)
!        res=dsigma_p_QL(p,rho,m,ga,A,r)*dp+dsigma_mu_QL(p,rho,m,ga,A,r)*dmu_r_QL(p,rho,m,ga,A,r)
    end function dsigma_r_QL

    !================================================================================
    ! Expansion coefficient of sigma near r = 0
    !================================================================================
    function sigma_2_QL(p,rho,ga,A) result(res)
        ! expansion coefficient at r = 0: sigma = sigma_0 + sigma_2 r^2 + ...
        use type_mod,only:pi
        real(wp),intent(in)::p,rho,ga,A
        real(wp)::res
        res = lam*p*rho*8*pi/3
    end function sigma_2_QL

end module anisotropy_QL_mod
