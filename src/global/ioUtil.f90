!io Utilities
!io sample: write(*,fmt='("<i=",i0,", reals=",2(g0,1x),">")')100,100.0,100.d0

module io_mod
    implicit none
    private
    public::fmt_es_use,fmt_title_use
    public::assert
    public::err_msg
    public::file_exists
    public::getUnit,writeTable,writeArray,reallocate_a,file_constTable,file_nCol,file_nRow,file_readTable
    public::to_char
    public::split_char
    public::save_bin,load_bin
    public::append_a
    public::reverse_a
    public::sort
    public::swap
    public::createFile
    public::arth
    public::alist
    public::linspace

    integer, parameter::dp=kind(1.d0)
    integer, parameter::wp=dp

    interface writeTable
        module procedure writeTable_R, writeTable_RR
    endinterface writeTable

    interface writeArray
        module procedure writeArrayHeader,writeTitle,writeArray_1,writeArray_2,writeArray_3,writeArray_4 &
        ,writeArray_5,writeArray_6,writeArray_7,writeArray_8,writeArray_9,writeArray_10,writeArray_11 &
        ,writeArray_12
    endinterface writeArray

    interface save_bin
        module procedure save_bin_v,save_bin_m,save_bin_m3,save_bin_m4
    end interface save_bin

    interface load_bin
        module procedure load_bin_v,load_bin_m,load_bin_m3,load_bin_m4
    end interface load_bin

    interface reallocate_a
        module procedure reallocate_iv_a,reallocate_rv_a, reallocate_rm_a,reallocate_cv_a, reallocate_cm_a
    endinterface reallocate_a

    interface append_a
        module procedure append_dv_a, append_cv_a
    end interface append_a

    interface reverse_a
        module procedure reverse_dv_a
    end interface reverse_a

    interface sort
        module procedure sort_dv_a
    end interface sort

    interface swap
        module procedure swap_d, masked_swap_ds
    end interface swap

    interface to_char
        module procedure wp_to_char, i_to_char
    end interface to_char

    interface arth
        module procedure arth_w
    end interface arth

    interface alist
        module procedure alist_wp
    end interface alist

    character(len=15),parameter::fmt_es_1='(15(es25.15))'
    character(len=15),parameter::fmt_es_2='(15(es16.8))'

    character(len=15),parameter::fmt_title_1='(15(a25))'
    character(len=15),parameter::fmt_title_2='(15(a16))'

    character(len=25)::fmt_es_use=trim(fmt_es_1)
    character(len=25)::fmt_title_use=trim(fmt_title_1)

contains
    function getUnit()
        implicit none
        integer:: getUnit
        integer:: i, ios
        logical:: lopen

        getUnit = 0
        do i = 1, 999
            inquire(unit = i, opened = lopen, iostat = ios)
            if ( ios == 0 ) then
                if ( .not. lopen ) then
                    getUnit = i
                    return
                end if
            end if
        end do
        write(*,*) 'getUnit: no free units'
        stop
    endfunction getUnit

    subroutine assert(cond,msg)
        logical,intent(in)::cond
        character(len=*),intent(in)::msg
        if (.not.cond) call err_msg(msg)
    end subroutine assert

    subroutine err_msg(msg)
        character(len=*)::msg
        write(*,*) msg
        stop
    end subroutine err_msg

    !Modified from NRUtil
    function reallocate_iv_a(p,n)
        integer, dimension(:), allocatable :: p, reallocate_iv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_iv_a(n),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_iv_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_iv_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_iv_a
    function reallocate_rv_a(p,n)
        real(wp), dimension(:), allocatable :: p, reallocate_rv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_rv_a(n),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_rv_a: problem in attempt to allocate memory')
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
        if (ierr /= 0) call err_msg('reallocate_rm_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm_a(1:min(nold,n),1:min(mold,m))=p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    endfunction reallocate_rm_a
    function reallocate_cv_a(p,n)
        complex(wp), dimension(:), allocatable :: p, reallocate_cv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_cv_a(n),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_cv_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_cv_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_cv_a
    function reallocate_cm_a(p,n,m)
        complex(wp), dimension(:,:), allocatable :: p, reallocate_cm_a
        integer, intent(in) :: n,m
        integer :: nold,mold,ierr
        allocate(reallocate_cm_a(n,m),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_cm_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_cm_a(1:min(nold,n),1:min(mold,m))=p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    endfunction reallocate_cm_a

    !Append array
    function append_dv_a(a_old,a_add)
        real(wp),dimension(:),intent(in)::a_old,a_add
        real(wp),dimension(:),allocatable::append_dv_a
        integer :: nold,nadd,ierr
        nold=size(a_old,1)
        nadd=size(a_add,1)
        allocate(append_dv_a(nold+nadd),stat=ierr)
        if (ierr /= 0) then
            write(*,*)'append_dv_a: problem in attempt to allocate memory'
            stop
        endif
        append_dv_a(1:nold)=a_old
        append_dv_a(nold+1:nold+nadd)=a_add
    endfunction append_dv_a
    function append_cv_a(a_old,a_add)
        complex(wp),dimension(:),intent(in)::a_old,a_add
        complex(wp),dimension(:),allocatable::append_cv_a
        integer :: nold,nadd,ierr
        nold=size(a_old,1)
        nadd=size(a_add,1)
        allocate(append_cv_a(nold+nadd),stat=ierr)
        if (ierr /= 0) then
            write(*,*)'append_cv_a: problem in attempt to allocate memory'
            stop
        endif
        append_cv_a(1:nold)=a_old
        append_cv_a(nold+1:nold+nadd)=a_add
    endfunction append_cv_a

    !reverse array
    subroutine reverse_dv_a(arr)
        real(wp),dimension(:),intent(inout)::arr
        real(wp),dimension(size(arr))::temp
        integer::i,n
        n=size(arr)
        do i=1,n
            temp(i)=arr(n+1-i)
        end do
        arr=temp
    end subroutine reverse_dv_a

    !sort array: From numerical recipes
    subroutine sort_dv_a(arr)
    real(wp), dimension(:), intent(inout) :: arr
    integer, parameter :: NN=15, NSTACK=50
    real(wp) :: a
    integer :: n,k,i,j,jstack,l,r
    integer, dimension(NSTACK) :: istack
    n=size(arr)
    jstack=0
    l=1
    r=n
    do
        if (r-l < NN) then
            do j=l+1,r
                a=arr(j)
                do i=j-1,l,-1
                    if (arr(i) <= a) exit
                    arr(i+1)=arr(i)
                end do
                arr(i+1)=a
            end do
            if (jstack == 0) return
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
        else
            k=(l+r)/2
            call swap(arr(k),arr(l+1))
            call swap(arr(l),arr(r),arr(l)>arr(r))
            call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
            call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
            i=l+1
            j=r
            a=arr(l+1)
            do
                do
                    i=i+1
                    if (arr(i) >= a) exit
                end do
                do
                    j=j-1
                    if (arr(j) <= a) exit
                end do
                if (j < i) exit
                call swap(arr(i),arr(j))
            end do
            arr(l+1)=arr(j)
            arr(j)=a
            jstack=jstack+2
            if (jstack > NSTACK) call err_msg('sort: NSTACK too small')
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
    end subroutine sort_dv_a

    !swap
    subroutine swap_d(a,b)
    real(wp), intent(inout) :: a,b
    real(wp) :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_d

    subroutine masked_swap_ds(a,b,mask)
    real(wp), intent(inout) :: a,b
    logical, intent(in) :: mask
    real(wp) :: swp
    if (mask) then
        swp=a
        a=b
        b=swp
    end if
    end subroutine masked_swap_ds

    !Read from data file
    subroutine file_nCol(path,nCol)
        implicit none
        character(len=*),intent(in)::path
        integer,intent(out)::nCol
        integer,parameter:: maxLineLen=1000,maxCol=30
        character(len=maxLineLen)::line, colTest(1:maxCol)
        integer::i, io, inFile
        inFile=getUnit()
        open(inFile,file=path,action='read',status='old')
        read(inFile,'(a)') line
        close(inFile)
        do i=1,maxCol
            read(line,*,iostat=io) colTest(1:i)
            if (io==-1) exit
        enddo
        nCol=i-1
    endsubroutine file_nCol

    subroutine file_nRow(path,nRow)
        implicit none
        character(len=*),intent(in)::path
        integer,intent(out)::nRow
        integer::io, inFile
        nRow=0
        inFile=getUnit()
        open(inFile,file=path,action='read',status='old')
        do
            read(inFile,*,iostat=io)
            if (io/=0) exit
            nRow=nRow+1
        enddo
        close(inFile)
    endsubroutine file_nRow

    subroutine file_readTable(path,table,skip) !read table with known nCol and nRow
        implicit none
        integer:: inFile
        character(len=*),intent(in):: path
        integer,optional,intent(in)::skip
        real(wp), intent(out):: table(1:,1:)
        integer:: i,nCol,nRow,io
        nCol=size(table,1)
        nRow=size(table,2)
        inFile=getUnit()
        open(inFile,file=path,status='old',action='read')
        if (present(skip)) then
            do i=1,skip
                read(inFile,*,iostat=io)
                if (io/=0) call err_msg('err: iostat - file_readTable')
            enddo
        endif
        do i=1,nRow
            read(inFile,*,iostat=io) table(1:nCol, i)
            if (io/=0) call err_msg('err: iostat - file_readTable')
        enddo
        close(inFile)
    endsubroutine file_readTable

    subroutine file_constTable(path,table,skip) !read and allocate table according to nCol and nRow
        implicit none
        character(len=*),intent(in)::path
        real(wp),allocatable,intent(inout)::table(:,:)
        integer,optional,intent(in)::skip
        integer::nCol,nRow
        if (allocated(table)) deallocate(table)
        call file_nCol(path,nCol)
        call file_nRow(path,nRow)
        if (present(skip))nRow=nRow-skip !title row & unit row
        allocate(table(1:nCol,1:nRow))
        call file_readTable(path,table,skip=skip)
    endsubroutine file_constTable

    !Writing data file
    subroutine writeTable_R(path, raw)
        implicit none
        character(len=*):: path
        real(wp), intent(in):: raw(1:, 1:)
        integer:: outFile
        integer:: nCol, nRow, i
        outFile=getUnit()
        open(outFile, file= path,status = 'replace', action='write')
        nCol=size(raw,1)
        nRow=size(raw,2)
        do i= 1, nRow
            write(outFile, fmt_es_use) raw(1:nCol,i)
        enddo
        close(outFile)
    endsubroutine writeTable_R

    subroutine writeTable_RR(path, raw1, raw2)
        implicit none
        character(len=*):: path
        real(wp), intent(in):: raw1(1:), raw2(1:, 1:)
        integer:: outFile
        integer:: nCol, nRow, i

        outFile=getUnit()
        open(outFile, file= path,status = 'replace', action='write')
        nCol=size(raw2,1)
        nRow=size(raw1,1)
        do i= 1, nRow
            write(outFile, fmt_es_use)raw1(i),  raw2(1:nCol,i)
        enddo
        close(outFile)
    endsubroutine writeTable_RR

    !Create file if not exist
    subroutine createFile(path,cc)
        character(len=*),intent(in)::path
        character(len=*),optional,intent(in):: cc
        integer::outFile
        logical::lexist
        integer:: i,pos,ncount

        outFile=getUnit()
        inquire(file=path,exist=lexist)
        if (.not. lexist) then
            open(outFile,file=path,status='replace')

            if (present(cc)) then
                ncount=0
                do i=1,len(cc)
                    if (cc(i:i)=='#') ncount=ncount+1
                end do
                if (ncount==0) then
                    write(outFile,*) cc
                    close(outFile)
                    return
                end if
                pos=0
                do i=1,ncount-1
                    pos=scan(cc(pos+1:),'#')+pos
                    write(outFile,'(a25)',advance='no') trim(cc(pos:scan(cc(pos+1:),'#')-1+pos))
                end do
                pos=scan(cc(pos+1:),'#')+pos
                write(outFile,'(a25)') trim(cc(pos:))
            end if

            close(outFile)
        endif
    end subroutine createFile

    function file_exists(path)
        character(len=*),intent(in)::path
        logical::file_exists
        inquire(file=path, exist=file_exists)
    end function file_exists

    !Write array
    subroutine writeArray_status(outFile,path,st,ac,po)
        integer,intent(in)::outFile
        character(len=*),intent(in)::path
        character(len=*),intent(in),optional::st,ac,po
        character(len=99)::status,access,position
        if (present(st)) then
            status=st
        else
            status='replace'
        end if
        if (present(ac)) then
            access=ac
        else
            access='sequential'
        end if
        if (present(po)) then
            position=po
        else
            position='asis'
        end if
        open(outFile, file= path,status=status, action='write',access=access,position=position)
    end subroutine writeArray_status

    subroutine writeArrayHeader(path,cc,st,ac,po)
        implicit none
        character(len=*),intent(in):: path,cc
        character(len=*),intent(in),optional::st,ac,po
        integer:: outFile,i,pos,ncount
        outFile=getUnit()
        call writeArray_status(outFile,path,st,ac,po)
        ncount=0
        do i=1,len(cc)
            if (cc(i:i)=='#') ncount=ncount+1
        end do
        if (ncount==0) then
            write(outFile,*) cc
            close(outFile)
            return
        end if
        pos=0
        do i=1,ncount-1
            pos=scan(cc(pos+1:),'#')+pos
            write(outFile,'(a25)',advance='no') trim(cc(pos:scan(cc(pos+1:),'#')-1+pos))
        end do
        pos=scan(cc(pos+1:),'#')+pos
        write(outFile,'(a25)') trim(cc(pos:))
        close(outFile)
    endsubroutine writeArrayHeader

    subroutine writeTitle(path,cc,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        character(len=*),dimension(:),intent(in):: cc
        character(len=*),intent(in),optional::st,ac,po
        integer:: outFile

        outFile=getUnit()
        call writeArray_status(outFile,path,st,ac,po)
        write(outFile,fmt_title_use) cc
        close(outFile)
    end subroutine writeTitle

    subroutine writeArray_1(path,a1,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i)
        enddo
        close(outFile)
    endsubroutine writeArray_1
    subroutine writeArray_2(path,a1,a2,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:), a2(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i)
        enddo
        close(outFile)
    endsubroutine writeArray_2
    subroutine writeArray_3(path,a1,a2,a3,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i)
        enddo
        close(outFile)
    endsubroutine writeArray_3
    subroutine writeArray_4(path,a1,a2,a3,a4,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i)
        enddo
        close(outFile)
    endsubroutine writeArray_4
    subroutine writeArray_5(path,a1,a2,a3,a4,a5,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i)
        enddo
        close(outFile)
    endsubroutine writeArray_5
    subroutine writeArray_6(path,a1,a2,a3,a4,a5,a6,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i)
        enddo
        close(outFile)
    endsubroutine writeArray_6
    subroutine writeArray_7(path,a1,a2,a3,a4,a5,a6,a7,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i)
        enddo
        close(outFile)
    endsubroutine writeArray_7
    subroutine writeArray_8(path,a1,a2,a3,a4,a5,a6,a7,a8,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i)
        enddo
        close(outFile)
    endsubroutine writeArray_8
    subroutine writeArray_9(path,a1,a2,a3,a4,a5,a6,a7,a8,a9,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i),a9(i)
        enddo
        close(outFile)
    endsubroutine writeArray_9
    subroutine writeArray_10(path,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:),a10(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i),a9(i),a10(i)
        enddo
        close(outFile)
    endsubroutine writeArray_10
    subroutine writeArray_11(path,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:),a10(:),a11(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i),a9(i),a10(i),a11(i)
        enddo
        close(outFile)
    endsubroutine writeArray_11
    subroutine writeArray_12(path,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,st,ac,po)
        implicit none
        character(len=*),intent(in):: path
        real(wp), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:),a10(:),a11(:),a12(:)
        character(len=*),intent(in),optional::st,ac,po
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac,po)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, fmt_es_use)a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i),a9(i),a10(i),a11(i),a12(i)
        enddo
        close(outFile)
    endsubroutine writeArray_12


    !===================================================================
    ! Saving data to binary file
    !===================================================================
    subroutine save_bin_v(filename,datav)
        character(len=*),intent(in)::filename
        real(wp),dimension(:),intent(in)::datav
        integer::uni

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='replace')
        write(uni) shape(datav)
        write(uni) datav
        close(uni)
    end subroutine save_bin_v

    subroutine save_bin_m(filename,datam)
        character(len=*),intent(in)::filename
        real(wp),dimension(:,:),intent(in)::datam
        integer::uni

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='replace')
        write(uni) shape(datam)
        write(uni) datam
        close(uni)
    end subroutine save_bin_m

    subroutine save_bin_m3(filename,datam)
        character(len=*),intent(in)::filename
        real(wp),dimension(:,:,:),intent(in)::datam
        integer::uni

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='replace')
        write(uni) shape(datam)
        write(uni) datam
        close(uni)
    end subroutine save_bin_m3

    subroutine save_bin_m4(filename,datam)
        character(len=*),intent(in)::filename
        real(wp),dimension(:,:,:,:),intent(in)::datam
        integer::uni

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='replace')
        write(uni) shape(datam)
        write(uni) datam
        close(uni)
    end subroutine save_bin_m4

    !===================================================================
    ! Loading data from binary file
    !===================================================================
    subroutine load_bin_v(filename,datav)
        character(len=*),intent(in)::filename
        real(wp),allocatable,dimension(:),intent(out)::datav
        integer::uni
        integer,dimension(1)::shape_v

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='old')
        read(uni) shape_v
        allocate(datav(shape_v(1)))
        read(uni) datav
        close(uni)
    end subroutine load_bin_v

    subroutine load_bin_m(filename,datam)
        character(len=*),intent(in)::filename
        real(wp),allocatable,dimension(:,:),intent(out)::datam
        integer::uni
        integer,dimension(2)::shape_m
        real(wp),allocatable,dimension(:)::tempv
!        integer::i,j

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='old')
        read(uni) shape_m
        allocate(tempv(shape_m(1)*shape_m(2)))
        read(uni) tempv
        datam = reshape(tempv,shape_m)
!        allocate(datam(shape_m(1),shape_m(2)))
!        do j=1,shape_m(2)
!            do i=1,shape_m(1)
!                read(uni) datam(i,j)
!            end do
!        end do
        close(uni)
    end subroutine load_bin_m

    subroutine load_bin_m3(filename,datam)
        character(len=*),intent(in)::filename
        real(wp),allocatable,dimension(:,:,:),intent(out)::datam
        integer::uni
        integer,dimension(3)::shape_m
        real(wp),allocatable,dimension(:)::tempv
!        integer::i,j,k

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='old')
        read(uni) shape_m
        allocate(tempv(shape_m(1)*shape_m(2)*shape_m(3)))
        read(uni) tempv
        datam = reshape(tempv,shape_m)
!        allocate(datam(shape_m(1),shape_m(2),shape_m(3)))
!        do k=1,shape_m(3)
!            do j=1,shape_m(2)
!                do i=1,shape_m(1)
!                    read(uni) datam(i,j,k)
!                end do
!            end do
!        end do
        close(uni)
    end subroutine load_bin_m3

    subroutine load_bin_m4(filename,datam)
        character(len=*),intent(in)::filename
        real(wp),allocatable,dimension(:,:,:,:),intent(out)::datam
        integer::uni
        integer,dimension(4)::shape_m
        real(wp),allocatable,dimension(:)::tempv
!        integer::i,j,k,l

        open(newunit=uni,file=filename,form='unformatted',access='stream',status='old')
        read(uni) shape_m
        allocate(tempv(shape_m(1)*shape_m(2)*shape_m(3)*shape_m(4)))
        read(uni) tempv
        datam = reshape(tempv,shape_m)
!        allocate(datam(shape_m(1),shape_m(2),shape_m(3),shape_m(4)))
!        do l=1,shape_m(4)
!            do k=1,shape_m(3)
!                do j=1,shape_m(2)
!                    do i=1,shape_m(1)
!                        read(uni) datam(i,j,k,l)
!                    end do
!                end do
!            end do
!        end do
        close(uni)
    end subroutine load_bin_m4

    !Convert to characters
    function i_to_char(i_in)
        integer,intent(in)::i_in
        character(len=:),allocatable::i_to_char
        allocate(character(len=25)::i_to_char)
        write(i_to_char,'(i0)')i_in
        i_to_char=trim(i_to_char)
    end function i_to_char

    function wp_to_char(wp_in,fmt_in,nospace)
        real(wp),intent(in)::wp_in
        character(len=*),intent(in),optional::fmt_in
        logical,intent(in),optional::nospace
        character(len=:),allocatable::fmt_use
        character(len=:),allocatable::wp_to_char
        allocate(character(len=25)::wp_to_char)
        if (present(fmt_in)) then
            fmt_use=fmt_in
        else
            fmt_use=fmt_es_use
        end if
        write(wp_to_char,fmt_use) wp_in
        wp_to_char=trim(wp_to_char)
        if (present(nospace)) then
            if (nospace) wp_to_char=trim(adjustl(wp_to_char))
        end if
    end function wp_to_char

    !Split string into array of characters
    function split_char(str,sep_char,le) result(char_array)
        character(len=*),intent(in)::str
        character(len=1),intent(in)::sep_char
        integer,intent(in)::le
        character(len=:),allocatable,dimension(:)::char_array
        integer,parameter::max_size=9999
        integer,dimension(max_size)::sep
        integer::nsep
        integer::i

        nsep=0
        do i=1,len(str) ! locate the positions of the sep_char
            if (str(i:i)==sep_char) then
                nsep=nsep+1
                sep(nsep)=i
            end if
        end do
        if (str(len(str):len(str)) /= sep_char) then ! padding the last character with sep_char
            nsep=nsep+1
            sep(nsep)=len(str)+1
        end if
        call assert(nsep>0,'err: str no separater: '//sep_char)
        allocate(character(len=le)::char_array(nsep))
        char_array(1)=str(1:sep(1)-1)
        do i=2,nsep
            char_array(i)=str(sep(i-1)+1:sep(i)-1)
        end do
    end function split_char

    !Gen array

    function arth_w(first,increment,n) !modified from numerical recipe nrutil
    real(wp), intent(in) :: first,increment
    integer, intent(in) :: n
    integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
    real(wp), dimension(n) :: arth_w
    integer :: k,k2
    real(wp) :: temp
    if (n > 0) arth_w(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth_w(k)=arth_w(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth_w(k)=arth_w(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth_w(k+1:min(k2,n))=temp+arth_w(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
    end function arth_w

    function alist_wp(ti,tf,nTot)
        real(wp),intent(in)::ti,tf
        integer,intent(in)::nTot
        real(wp),allocatable::alist_wp(:)
        integer::i
        allocate(alist_wp(nTot))
        do i=1,nTot
            alist_wp(i)=ti+i*(tf-ti)/nTot
        end do
    end function alist_wp

    function linspace(from, to, n) result(array)
        real(wp), intent(in) :: from, to
        integer:: n
        real(wp), allocatable, dimension(:) :: array
        real(wp) :: range
        integer :: i

        allocate(array(n))
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if

        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end function linspace

endmodule io_mod
