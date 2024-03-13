module anisotropy_type_mod
    use type_mod,only:wp
    implicit none
    private
    public:: set_par, sigma_func, dsigma_func, sigma_coef
    public:: aniso_type_var
!    public:: assign_aniso_type

    interface
        subroutine set_par(x)
            import wp
            real(wp),intent(in)::x
        end subroutine set_par

        function sigma_func(p,rho,m,r) result(res)
            import wp
            real(wp),intent(in)::p,rho,m,r
            real(wp)::res
        end function sigma_func

        function dsigma_func(p,rho,m,ga,A,r) result(res)
            import wp
            real(wp),intent(in)::p,rho,m,ga,A,r
            real(wp)::res
        end function dsigma_func

        function sigma_coef(p,rho,ga,A) result(res)
            import wp
            real(wp),intent(in)::p,rho,ga,A
            real(wp)::res
        end function sigma_coef

    end interface

    type aniso_type_var
        procedure(set_par), pointer, nopass :: set

        procedure(sigma_func), pointer, nopass :: s
        procedure(dsigma_func), pointer, nopass :: dsp
        procedure(dsigma_func), pointer, nopass :: dsrho
        procedure(dsigma_func), pointer, nopass :: dsmu

        procedure(dsigma_func), pointer, nopass :: ddsp
        procedure(dsigma_func), pointer, nopass :: ddspmu
        procedure(dsigma_func), pointer, nopass :: dmu
        procedure(dsigma_func), pointer, nopass :: ds

        procedure(sigma_coef), pointer, nopass :: s2
    end type aniso_type_var

!    contains
!
!    subroutine assign_aniso_type(ain,aout)
!        type(aniso_type_var),intent(in)::ain
!        type(aniso_type_var),intent(out)::aout
!
!        aout%set        =>       ain%set
!
!        aout%s          =>       ain%s
!        aout%dsp        =>       ain%dsp
!        aout%dsmu       =>       ain%dsmu
!
!        aout%ddsp       =>       ain%ddsp
!        aout%ddspmu     =>       ain%ddspmu
!        aout%dmu        =>       ain%dmu
!
!        aout%s2         =>       ain%s2
!    end subroutine assign_aniso_type

end module anisotropy_type_mod
