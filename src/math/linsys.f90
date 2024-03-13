module linsys_mod
    implicit none
    private
    public::det,tri_matrix,LinSys_Cramer
    interface det
        module procedure det_r, det_c
    end interface
    interface tri_matrix
        module procedure tri_matrix_r, tri_matrix_c
    end interface
    interface LinSys_Cramer
        module procedure LinSys_Cramer_r, LinSys_Cramer_c
    end interface

contains

    function det_r(a)
        real(8),intent(in)::a(:,:)
        real(8) :: det_r, b(size(a,1),size(a,2)), output
        integer :: i,l
        if (size(a,1)/=size(a,2)) call err_msg('err: a is not square matrix')
        call tri_matrix(a,b,l)
        output = l                      !   l represents the sign change due to swapping 2 rows;
        do i= 1,size(a,1)
            output = output*b(i,i)
        enddo
        det_r=output
    end function

    function det_c(a)
        implicit none
        complex(8),intent(in)::a(:,:)
        complex(8) :: det_c, b(size(a,1),size(a,2)), output
        integer :: i,l
        if (size(a,1)/=size(a,2)) call err_msg('err: a is not square matrix')
        call tri_matrix(a,b,l)
        output = l                      !   l represents the sign change due to swapping 2 rows;
        do i= 1,size(a,1)
            output = output*b(i,i)
        enddo
        det_c=output
    end function

    !************************************************************************************
    !   converting a matrix to upper triangular form
    !************************************************************************************
    subroutine tri_matrix_r(a, b, l)
        implicit none
        real(8) :: a(:, :), b(:,:)
        real(8) :: factor
        integer :: i, j, n, k, l
        logical :: det_exist

        det_exist = .true.
        !!parti : initialization
        n=size(a,1)
        l = 1
        do i = 1,n
            do j = 1,n
                b(i,j) = a(i,j)
            enddo
        enddo
        !!partii :
        do i = 1,n-1
            !!partiia : check existence of tri_matrix
            if (b(i,i) == 0.d0) then
                det_exist = .false.
                do j = i+1, n
                    if (b(j,i) /= 0.d0) then
                        det_exist = .true.
                        do k = 1, n
                            call swap(b(i,k), b(j,k))
                        enddo
                        l = -l
                        exit
                    endif
                enddo
            endif
            if (.not. det_exist) call err_msg("err: tri_matrix; the matrix cannot be changed to ut form")
            !!partiib : make triangular matrix
            do j = i+1, n
                factor = b(j,i)/b(i,i)
                do k = i, n
                    b(j,k) = b(j,k) - factor * b(i,k)
                enddo
            enddo
        enddo
    contains
        subroutine swap(a,b)
            real(8) :: a, b, c
            c=a
            a=b
            b=c
        endsubroutine swap
    endsubroutine tri_matrix_r

    subroutine tri_matrix_c(a, b, l)
        implicit none
        complex(8) :: a(:, :), b(:,:)
        complex(8) :: factor
        integer :: i, j, n, k, l
        logical :: det_exist

        det_exist = .true.
        !!parti : initialization
        n=size(a,1)
        l = 1
        do i = 1,n
            do j = 1,n
                b(i,j) = a(i,j)
            enddo
        enddo
        !!partii :
        do i = 1,n-1
            !!partiia : check existence of tri_matrix
            if (b(i,i) == 0.d0) then
                det_exist = .false.
                do j = i+1, n
                    if (b(j,i) /= 0.d0) then
                        det_exist = .true.
                        do k = 1, n
                            call swap(b(i,k), b(j,k))
                        enddo
                        l = -l
                        exit
                    endif
                enddo
            endif
            if (.not. det_exist) call err_msg("err: tri_matrix; the matrix cannot be changed to ut form")
            !!partiib : make triangular matrix
            do j = i+1, n
                factor = b(j,i)/b(i,i)
                do k = i, n
                    b(j,k) = b(j,k) - factor * b(i,k)
                enddo
            enddo
        enddo
    contains
        subroutine swap(a,b)
            complex(8) :: a, b, c
            c=a
            a=b
            b=c
        endsubroutine swap
    endsubroutine tri_matrix_c

    !************************************************************************************
    !   Solving Linear System: Cramer's Rule
    !************************************************************************************
    !   Works quite well for a 4x4 system
    !   Gives huge error for a 6s6 system (5~10%)
    subroutine LinSys_Cramer_r(a,b,vec)
        !   Solving system "[a] vec = b"
        implicit none
        real(8),intent(in) :: a(:,:),b(:)
        real(8),intent(out):: vec(:)
        real(8) :: M(size(a,1),size(a,1))
        integer :: i
        if (det(a) == 0.d0) call err_msg("err: LinSys_Cramer; det = 0")
        do i = 1,size(a,1)
            M=a
            M(:,i) = b(:)
            vec(i)=det(M)/det(a)
        enddo
    endsubroutine LinSys_Cramer_r
    subroutine LinSys_Cramer_c(a, b,vec)
        !   Solving system "[a] vec = b"
        implicit none
        complex(8),intent(in) :: a(:,:),b(:)
        complex(8),intent(out):: vec(:)
        complex(8)::M(size(a,1),size(a,1))
        integer :: i
        if (det(a) == 0.d0) call err_msg("err: LinSys_Cramer; det = 0")
        do i = 1,size(a,1)
            M=a
            M(:,i) = b(:)
            vec(i)=det(M)/det(a)
        enddo
    endsubroutine LinSys_Cramer_c

    subroutine err_msg(msg)
        character(len=*)::msg
        write(*,*) msg
        stop
    end subroutine err_msg
endmodule linsys_mod
