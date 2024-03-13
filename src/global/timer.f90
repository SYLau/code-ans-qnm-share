module ctime_mod
    implicit none
    private
    public::ctimer
    public::print_time
    real(8)::t1,t2
    contains
    subroutine ctimer(n,msg)
        integer,intent(in)::n
        character(len=*),intent(in),optional::msg
        if (present(msg)) print*,">>>Timing: ",msg
        if (n==1) call cpu_time(t1)
        if (n==2) then
            call cpu_time(t2)
            print("(1a, f15.8)"), "Elapsed CPU time (s) = ", t2-t1
            !print("(x,1a, f25.16)"), "Elapsed CPU time (s) = ", t2-t1
        endif
    endsubroutine ctimer

    subroutine print_time(msg)
        character(len=*),intent(in),optional::msg
        real(8)::t
        call cpu_time(t)
        if (present(msg)) then
            print("(1a, f15.8, 1a)"),">>> "//trim(msg)//" [Time (s) = ",t," ]"
        else
            print("(1a, f15.8, 1a)"),">>> [Time (s) = ",t," ]"
        endif
    end subroutine print_time
endmodule ctime_mod
