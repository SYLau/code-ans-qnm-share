program Sid_LD
    use type_mod,only:wp
    use dirUtil_mod,only:dir_create
    implicit none

    real(wp)::rho0
    real(wp)::lam
    character(len=:),allocatable::eos_file
    real(wp)::epsmin_tov
    real(wp)::ell
    character(len=:),allocatable::aniso_EOS
    namelist /ld_single_par/ rho0,lam,eos_file,epsmin_tov,ell,aniso_EOS

    integer::n
    real(wp)::f1,f2
    integer::display_Ain

    real(wp)::fg,fwidth
    real(wp)::fi
    character(len=:),allocatable::save_mode

    namelist /mode_search_ctrl/ n,f1,f2,display_Ain,fg,fwidth,save_mode,fi
    real(wp),allocatable,dimension(:)::out_f,out_A
    character(len=:),allocatable::option

    !======================================================================================================================================
    ! Create output folder
    !======================================================================================================================================
    call dir_create('results_single')

    !======================================================================================================================================
    ! Initialize global arrays
    !======================================================================================================================================
    allocate(out_f(0),out_A(0))

    !======================================================================================================================================
    initialize_output_file: block
        integer::uni1,uni2

        open(newunit=uni1,file='results_single/freq.dat',status='replace',action='write')
        write(uni1,'(7a16)') [character(len=12)::'mode', 'Re(f) (Hz)', 'Im(f) (Hz)', 'fguess (Hz)', 'fwidth', 'rho0', 'lam']
        open(newunit=uni2,file='results_single/mode.dat',status='replace',action='write')
        close(uni1)
        close(uni2)
    end block initialize_output_file


    !======================================================================================================================================
    main_loop:do

        !======================================================================================================================================
        choose_options: block
            use ctime_mod,only: print_time

            write(*,*) '-------------------------------------'
            call print_time()
            if (allocated(option)) deallocate(option)
            allocate(character(len=10)::option)
            write(*,*) '-------------------------------------'
            write(*,*) 'Change parameters in inlist_ld_single.nml'
            write(*,*) 'Choose from the following options'
            write(*,*) '1. Scan A_in'
            write(*,*) '2. Solve for quasi-normal mode'
            write(*,*) '0. quit'
            write(*,*) '-------------------------------------'
            read(*,*) option
            option=trim(option)
            if (option == '0') exit
        end block choose_options

        !======================================================================================================================================
        read_inlist: block
            integer:: uni

            if (allocated(eos_file)) deallocate(eos_file)
            if (allocated(aniso_EOS)) deallocate(aniso_EOS)
            if (allocated(save_mode)) deallocate(save_mode)
            allocate(character(len=30)::eos_file)
            allocate(character(len=30)::aniso_EOS)
            allocate(character(len=30)::save_mode)
            open(newunit=uni,file='inlist_ld_single.nml',status='old')
            read(uni,nml=ld_single_par)
            read(uni,nml=mode_search_ctrl)
            close(uni)
            eos_file=trim(eos_file)
            aniso_EOS=trim(aniso_EOS)
            save_mode=trim(save_mode)
        end block read_inlist

        !======================================================================================================================================
        solve_problem: block
            use type_mod,only:M_sun
            use eos_mod,only:eos_pt,eos_read,eos_p,eos_rho,eos_ga_r_sm
            use tov_aniso_pt_mod,only:tov_aniso_pt_bg
            use tov_aniso_pt_mod,only:solve_tov_aniso_pt
            use anisotropy_type_mod,only:aniso_type_var
            use bg_mod,only:load_bgs
            use io_mod,only:writeArray
            use io_mod,only:fmt_title_use,fmt_es_use

            type(tov_aniso_pt_bg)::tovs
            real(wp),allocatable,dimension(:)::pt

            type(aniso_type_var)::aEOS      !Anisotropic equation of state

            !======================================================================================================================================
            assign_EOS: block
                use anisotropy_QL_mod,only:aniso_QL
                use anisotropy_QL_2_mod,only:aniso_QL2
                use anisotropy_BL_mod,only:aniso_BL

                type(aniso_QL)::aQL
                type(aniso_QL2)::aQL2
                type(aniso_BL)::aBL

                if (associated(aEOS%set)) nullify(aEOS%set,aEOS%s,aEOS%dsp,aEOS%dsrho,aEOS%dsmu,aEOS%dmu,aEOS%s2)
                select case(aniso_EOS)
                    case ('QL')
                        aEOS%set        =>       aQL%set
                        aEOS%s          =>       aQL%s
                        aEOS%dsp        =>       aQL%dsp
                        aEOS%dsrho      =>       aQL%dsrho
                        aEOS%dsmu       =>       aQL%dsmu
                        aEOS%dmu        =>       aQL%dmu
                        aEOS%ds         =>       aQL%ds
                        aEOS%s2         =>       aQL%s2
                    case ('QL2')
                        aEOS%set        =>       aQL2%set
                        aEOS%s          =>       aQL2%s
                        aEOS%dsp        =>       aQL2%dsp
                        aEOS%dsrho      =>       aQL2%dsrho
                        aEOS%dsmu       =>       aQL2%dsmu
                        aEOS%dmu        =>       aQL2%dmu
                        aEOS%ds         =>       aQL2%ds
                        aEOS%s2         =>       aQL2%s2
                    case ('BL')
                        aEOS%set        =>       aBL%set
                        aEOS%s          =>       aBL%s
                        aEOS%dsp        =>       aBL%dsp
                        aEOS%dsrho      =>       aBL%dsrho
                        aEOS%dsmu       =>       aBL%dsmu
                        aEOS%dmu        =>       aBL%dmu
                        aEOS%ds         =>       aBL%ds
                        aEOS%s2         =>       aBL%s2
                    case default
                        error stop
                end select

            end block assign_EOS

            !Set format for writeArray
            fmt_es_use= '(10(es16.8))'
            fmt_title_use='(10(a16))'
            !Set the value of Lambda in aniso_QL
            call aEOS%set(lam)

            call eos_read(eos_file)  !Read the EOS table
            pt=eos_pt(rho0)
            call solve_tov_aniso_pt(eos_p(rho0),pt,epsmin_tov,eos_rho,eos_ga_r_sm,aEOS%s,tovs)

            output_tov: block
                use type_mod, only: G,c,pi
                real(wp),allocatable,dimension(:)::cs2,ct2
                real(wp)::dp,eLam
                integer::i

                cs2 = tovs%ga*tovs%p/(tovs%rho*c**2+tovs%p)
                if (allocated(ct2)) deallocate(ct2)
                allocate(ct2(size(cs2)))
                do i=1,size(cs2)
                    eLam = 1._wp/(1-2*G*tovs%m(i)/c**2/tovs%r(i))
                    dp = -G/c**2/tovs%r(i)**2*(tovs%rho(i)+tovs%p(i)/c**2) &
                    *(G/c**2*tovs%m(i)+G/c**4*pi*4*tovs%r(i)**3*tovs%p(i))*eLam
                    ct2(i) = cs2(i)*(1- aEOS%ds(G/c**4*tovs%p(i),G/c**2*tovs%rho(i),G/c**2*tovs%m(i) &
                    ,tovs%ga(i),tovs%A(i),tovs%r(i))/dp )
                end do

                call writeArray('results_single/tov_soln.txt',[character(len=12)::'r (cm)', 'P (cgs)', 'rho (cgs)' &
                , 'm (g)', 'nu', 'ga', 's (cm^{-2})', 'cs^2', 'ct^2'])
                call writeArray('results_single/tov_soln.txt',tovs%r,tovs%p,tovs%rho,tovs%m,tovs%nu,tovs%ga,tovs%s &
                , cs2, ct2,st='old',po='append')
            end block output_tov

!            call writeArray('results_single/tov_soln.txt',[character(len=12)::'r (cm)', 'P (cgs)', 'rho (cgs)' &
!            , 'm (g)', 'nu', 'ga', 's (cm^{-2})'])
!            call writeArray('results_single/tov_soln.txt',tovs%r,tovs%p,tovs%rho,tovs%m,tovs%nu,tovs%ga,tovs%s &
!            ,st='old',po='append')

            print*, 'M = ', tovs%m(tovs%n)/M_sun

            call load_bgs(tovs%r,tovs%p,tovs%rho,tovs%m,tovs%nu,tovs%ga,tovs%A,tovs%pt_pos)

            !======================================================================================================================================
            set_r_initial: block
                use puls_ld_aniso_mod,only:n_initial
                integer::i
                real(wp),parameter::initial=1.e-2

                do i=1,tovs%n-1
                    if (tovs%r(i)<=tovs%r(tovs%n)*initial .and. tovs%r(i+1)>tovs%r(tovs%n)*initial) exit
                end do
                n_initial = i
            end block set_r_initial

            !======================================================================================================================================
            solve_full_gr: block
                use eigen_ld_aniso_mod,only:scan_Ain_r,get_ld_var,pul_sol
                use eigen_ld_aniso_mod,only:get_eigenmodes_ld
                use ld_eq_aniso_mod,only:def_sigma
                use ctime_mod,only:print_time

                complex(wp):: fguess,fguess1,fguess2

                complex(wp),allocatable,dimension(:)::plt_f,plt_y
                complex(wp),allocatable,dimension(:)::sol_f
                type(pul_sol),allocatable,dimension(:)::sol_puls
                integer::ostat

                call def_sigma(aEOS%s,aEOS%dsp,aEOS%dsrho,aEOS%dsmu,aEOS%ds,aEOS%s2)

                if (option == '1') then
                    call print_time(msg='scanning A_{in}')

                    call scan_Ain_r(ell,f1,f2,n,plt_f,plt_y, fi=fi) !Calculate amplitude of ingoing wave as function of Re(f) (Hz) to estimate the position of the mode

                    !============================================================================================================
                    output_Ain: block
                        use io_mod,only:writeArray
                        use io_mod,only:append_a
                        use sort_nr_mod,only:sort
                        integer::uni
                        integer::i

                        call writeArray('results_single/ingoing_amp.txt',[character(len=12)::'f (Hz)', 'A_in'])
                        call writeArray('results_single/ingoing_amp.txt',real(plt_f),abs(plt_y),st='old',po='append')

                        !===========================================================================================================
                        plot_results:block
                            use dirUtil_mod,only:check_type

                            select case (display_Ain)
                            case(1)
                                if (check_type()=='cmd') &
                                call execute_command_line('cd results_single && start wgnuplot -p plot_amp.gp')
                                if (check_type()=='bash') &
                                call execute_command_line('cd results_single && gnuplot -p plot_amp.gp &')
                            case(2)
                                if (check_type()=='cmd') &
                                call execute_command_line('cd results_single && start wgnuplot -p plot_amp_full.gp')
                                if (check_type()=='bash') &
                                call execute_command_line('cd results_single && gnuplot -p plot_amp_full.gp &')
                            end select
!                            if (display_Ain==1) then
!                                if (check_type()=='cmd') &
!                                call execute_command_line('cd results_single && start wgnuplot -p plot_amp.gp')
!                                if (check_type()=='bash') &
!                                call execute_command_line('cd results_single && gnuplot -p plot_amp.gp &')
!                            end if
                        end block plot_results

                        out_f = append_a(out_f,real(plt_f))
                        out_A = append_a(out_A,abs(plt_y))

                        call sort(out_f, out_A)

                        open(newunit=uni, file = 'results_single/ingoing_amp_full.txt', status='replace', action = 'write')
                        write(uni, '(2a16)') 'f (Hz)', 'A_in'
                        do i=1,size(out_f)
                            write(uni, '(2es16.8)') out_f(i), out_A(i)
                        end do
                        close(uni)
                    end block output_Ain

                elseif (option == '2') then

                    if (allocated(sol_f)) deallocate(sol_f)
                    if (allocated(sol_puls)) deallocate(sol_puls)
                    allocate(sol_f(1),sol_puls(1))


                    call print_time(msg='solving for eigenmode')
                    fguess=cmplx(fg,fi,kind=wp)
                    fguess1=fguess*(1+fwidth)
                    fguess2=fguess*(1-fwidth)
                    print*, 'fguess = ', fg
                    print*, 'fwidth = ', fwidth
                    call get_eigenmodes_ld(ell,fguess,fguess1,fguess2,sol_f(1),sol_puls(1) &
                    ,quit_if_error=.false.,ostat=ostat)
                    if (ostat /= 0) exit solve_full_gr

                    call get_ld_var(ell,sol_f,sol_puls) ! convert the solutions in H1, K, W, X, H0, V and store in puls

                    print*, 'eigen-frequency (kHz) = ', real(sol_f)/1.e3
                    print*, 'imaginary part (Hz) = ', imag(sol_f)

                    !======================================================================================================================================
                    output_results: block
                        integer::i
                        real(wp),allocatable,dimension(:)::plt1,plt2,plt3,plt4,plt5,plt6,plt7!,plt8
                        integer::uni1, uni2

                        open(newunit=uni1,file='results_single/freq.dat',status='old',action='write',position='append')
                        open(newunit=uni2,file='results_single/mode.dat',status='old',action='write',position='append')

                        write(uni1,'(a16,6es16.8)') save_mode, real(sol_f(1)),imag(sol_f(1)), fg, fwidth, rho0, lam
                        write(uni2,'(a)') '################################################'
                        write(uni2,'(a)') '"'//save_mode//'"'
    !                    write(uni2,'(9a16)')'r (km)', 'Re(H1)', 'Re(K)', 'Re(W)', 'Re(X)', 'Im(H1)', 'Im(K)', 'Im(W)', 'Im(X)'
                        write(uni2,'(8a16)')[character(len=12)::'r (km)' &
                        , 'Re(H1)', 'Re(K)', 'Re(W)', 'Re(X)', 'Re(H0)', 'Re(V)', 'Re(Z)']

                        plt1=real(sol_puls(1)%H1); plt2=real(sol_puls(1)%K)
                        plt3=real(sol_puls(1)%W);  plt4=real(sol_puls(1)%Xp)
                        plt5=real(sol_puls(1)%H0); plt6=real(sol_puls(1)%V)
                        plt7=real(sol_puls(1)%Z)

                        do i=1,size(sol_puls(1)%t),6
                            write(uni2,'(8es16.8)') sol_puls(1)%t(i), plt1(i),plt2(i),plt3(i),plt4(i) &
                            ,plt5(i),plt6(i),plt7(i)
                        end do
                        write(uni2,'(a)') new_line('a')

                        close(uni1)
                        close(uni2)
                    end block output_results
                else
                    print*, 'invalid input'
                    exit solve_full_gr
                end if

            end block solve_full_gr

        end block solve_problem

    end do main_loop

end program Sid_LD
