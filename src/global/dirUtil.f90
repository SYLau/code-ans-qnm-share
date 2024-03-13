module dirUtil_mod
    implicit none
    private
    public::check_type
    public::dir_create
    public::dir_clear
    public::dir_remove

    contains

    function check_type() result(term)
        character(len=:),allocatable::term
        integer::env_stat

            call get_environment_variable('ComSpec',status=env_stat)

            if (env_stat==0) then
                term = 'cmd'
            elseif (env_stat==1) then
                term = 'bash'
            else
                error stop 'err: check_type unidentified error'
            end if

    end function check_type

    subroutine dir_create(dir)
        character(len=*),intent(in)::dir
        character(len=:),allocatable::term

        term = check_type()

        select case(term)
            case ('cmd')
                call execute_command_line('if not exist ".\'//trim(dir)//'" mkdir "'//trim(dir)//'"') !create new directory
            case ('bash')
                call execute_command_line('mkdir -p "'//trim(dir)//'"')
            case default
                print*, term
                stop
        end select


    end subroutine dir_create

    subroutine dir_clear(dir)
        character(len=*)::dir
        character(len=:),allocatable::term

        term = check_type()

        select case(term)
            case ('cmd')
                call execute_command_line('cd "'//trim(dir)//'" && DEL /F/Q/S *.* 2>NUL') !delete all files in output
                call execute_command_line('cd ..')
            case ('bash')
                call execute_command_line('cd "'//trim(dir)//'" && rm -f *.*')
                call execute_command_line('cd ..')
            case default
                print*, term
                stop
        end select

    end subroutine dir_clear

    subroutine dir_remove(dir)
        character(len=*)::dir
        character(len=:),allocatable::term

        term = check_type()

        select case(term)
            case ('cmd')
                call execute_command_line('rmdir /Q/S "'//trim(dir)//'" 2>Nul') !remove directory and all files within
            case ('bash')
                call execute_command_line('rm -rf "'//trim(dir)//'"') !-r: remove the directory; -f force, suppress warnings
            case default
                print*, term
                stop
        end select

    end subroutine dir_remove

end module dirUtil_mod
