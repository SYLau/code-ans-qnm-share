module sort_nr_mod
    use type_mod,only:wp
    implicit none
    private
    public::sort
    public::indexx

    interface sort
        module procedure sort1_wp,sort2_wp,sort3_wp,sort1_i,sort2_i,sort3_i
    end interface sort

    interface indexx
        module procedure indexx_wp, indexx_i
    end interface indexx

    contains

    !===============================================================================================
    ! sorting real(wp)
    !===============================================================================================
    subroutine sort1_wp(arr)
        implicit none
        real(wp), dimension(:), intent(inout) :: arr
        integer :: ndum
        integer, allocatable, dimension(:) :: index
        call indexx(arr,index)
        arr=arr(index)
    end subroutine sort1_wp

    subroutine sort2_wp(arr,slave1)
        implicit none
        real(wp), dimension(:), intent(inout) :: arr,slave1
        integer :: ndum
        integer, allocatable, dimension(:) :: index
        if (size(arr)/=size(slave1)) error stop 'sort2_wp'
        call indexx(arr,index)
        arr=arr(index)
        slave1=slave1(index)
    end subroutine sort2_wp

    subroutine sort3_wp(arr,slave1,slave2)
        implicit none
        real(wp), dimension(:), intent(inout) :: arr,slave1,slave2
        integer :: ndum
        integer, allocatable, dimension(:) :: index
        if (size(arr)/=size(slave1) .or. size(arr)/=size(slave2)) error stop 'sort3_wp'
        call indexx(arr,index)
        arr=arr(index)
        slave1=slave1(index)
        slave2=slave2(index)
    end subroutine sort3_wp

    !===============================================================================================
    ! sorting integer
    !===============================================================================================

    subroutine sort1_i(arr)
        implicit none
        integer, dimension(:), intent(inout) :: arr
        integer :: ndum
        integer, allocatable, dimension(:) :: index
        call indexx(arr,index)
        arr=arr(index)
    end subroutine sort1_i

    subroutine sort2_i(arr,slave1)
        implicit none
        integer, dimension(:), intent(inout) :: arr,slave1
        integer :: ndum
        integer, allocatable, dimension(:) :: index
        if (size(arr)/=size(slave1)) error stop 'sort2_i'
        call indexx(arr,index)
        arr=arr(index)
        slave1=slave1(index)
    end subroutine sort2_i

    subroutine sort3_i(arr,slave1,slave2)
        implicit none
        integer, dimension(:), intent(inout) :: arr,slave1,slave2
        integer :: ndum
        integer, allocatable, dimension(:) :: index
        if (size(arr)/=size(slave1) .or. size(arr)/=size(slave2)) error stop 'sort3_i'
        call indexx(arr,index)
        arr=arr(index)
        slave1=slave1(index)
        slave2=slave2(index)
    end subroutine sort3_i


    !===============================================================================================
    ! Indexx routine for real (wp)
    !===============================================================================================
    subroutine indexx_wp(arr,index)
    implicit none
    real(wp), dimension(:), intent(in) :: arr
    integer, allocatable, dimension(:), intent(out) :: index
    integer, parameter :: NN=15, NSTACK=50
    real(wp) :: a
    integer :: n,k,i,j,indext,jstack,l,r
    integer, dimension(NSTACK) :: istack

    n=size(arr)
    allocate(index(n))
    index=arth_i(1,1,n)
    jstack=0
    l=1
    r=n
    do
        if (r-l < NN) then
            do j=l+1,r
                indext=index(j)
                a=arr(indext)
                do i=j-1,l,-1
                    if (arr(index(i)) <= a) exit
                    index(i+1)=index(i)
                end do
                index(i+1)=indext
            end do
            if (jstack == 0) return
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
        else
            k=(l+r)/2
            call swap_i(index(k),index(l+1))
            call icomp_xchg(index(l),index(r))
            call icomp_xchg(index(l+1),index(r))
            call icomp_xchg(index(l),index(l+1))
            i=l+1
            j=r
            indext=index(l+1)
            a=arr(indext)
            do
                do
                    i=i+1
                    if (arr(index(i)) >= a) exit
                end do
                do
                    j=j-1
                    if (arr(index(j)) <= a) exit
                end do
                if (j < i) exit
                call swap_i(index(i),index(j))
            end do
            index(l+1)=index(j)
            index(j)=indext
            jstack=jstack+2
            if (jstack > NSTACK) error stop 'indexx: NSTACK too small'
            if (r-i+1 >= j-l) then
                istack(jstack)=r
                istack(jstack-1)=i
                r=j-1
            else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
            end if
        end if
    end do
    contains
!BL
        subroutine icomp_xchg(i,j)
        integer, intent(inout) :: i,j
        integer :: swp
        if (arr(j) < arr(i)) then
            swp=i
            i=j
            j=swp
        end if
        end subroutine icomp_xchg
    end subroutine indexx_wp

    !===============================================================================================
    ! indexx routine for integers
    !===============================================================================================

    subroutine indexx_i(iarr,index)
    implicit none
    integer, dimension(:), intent(in) :: iarr
    integer, allocatable, dimension(:), intent(out) :: index
    integer, parameter :: NN=15, NSTACK=50
    integer :: a
    integer :: n,k,i,j,indext,jstack,l,r
    integer, dimension(NSTACK) :: istack

    n=size(iarr)
    allocate(index(n))
    index=arth_i(1,1,n)
    jstack=0
    l=1
    r=n
    do
        if (r-l < NN) then
            do j=l+1,r
                indext=index(j)
                a=iarr(indext)
                do i=j-1,l,-1
                    if (iarr(index(i)) <= a) exit
                    index(i+1)=index(i)
                end do
                index(i+1)=indext
            end do
            if (jstack == 0) return
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
        else
            k=(l+r)/2
            call swap_i(index(k),index(l+1))
            call icomp_xchg(index(l),index(r))
            call icomp_xchg(index(l+1),index(r))
            call icomp_xchg(index(l),index(l+1))
            i=l+1
            j=r
            indext=index(l+1)
            a=iarr(indext)
            do
                do
                    i=i+1
                    if (iarr(index(i)) >= a) exit
                end do
                do
                    j=j-1
                    if (iarr(index(j)) <= a) exit
                end do
                if (j < i) exit
                call swap_i(index(i),index(j))
            end do
            index(l+1)=index(j)
            index(j)=indext
            jstack=jstack+2
            if (jstack > NSTACK) error stop 'indexx: NSTACK too small'
            if (r-i+1 >= j-l) then
                istack(jstack)=r
                istack(jstack-1)=i
                r=j-1
            else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
            end if
        end if
    end do

    contains
        subroutine icomp_xchg(i,j)
        integer, intent(inout) :: i,j
        integer :: swp
        if (iarr(j) < iarr(i)) then
            swp=i
            i=j
            j=swp
        end if
        end subroutine icomp_xchg
    end subroutine indexx_i

    !===============================================================================================
    ! arth routine
    !===============================================================================================
    function arth_i(first,increment,n)
    integer, intent(in) :: first,increment,n
    integer, dimension(n) :: arth_i
    integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
    integer :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth_i(k)=arth_i(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
    end function arth_i

    !===============================================================================================
    ! swap routine
    !===============================================================================================

    subroutine swap_i(a,b)
    integer, intent(inout) :: a,b
    integer :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_i

end module sort_nr_mod
