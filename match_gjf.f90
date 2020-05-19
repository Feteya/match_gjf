!! Input : 1.gjf , 2.gjf match_input
!! 将2.gjf的原子序号重新编号，让两个结构的原子序号一致
program match_gjf
    implicit none

    type :: atoms
        integer             :: natoms           !! 原子总个数
        integer             :: charge           !! 电荷
        integer             :: spin_mult        !! 自旋多重度
        integer,allocatable :: symbols(:)       !! 元素符号
        real,allocatable    :: coordinate(:,:)  !! 坐标 3*natoms
        real,allocatable    :: connect(:,:)     !! 链接度 n*n
        integer,allocatable :: link(:)          !! 共链接几个原子
    end type atoms

    !! 命令行参数组
    integer                         :: narg
    character(len=256),allocatable  :: arg(:)
    
    !! 两个分子的分子结构
    type(atoms)                     :: atoms1, atoms2

    !! 输出链接度信息时的临时变量
    character(len=120)              :: conn
    character(len=120)              :: temp_string

    !! 重新编号所用到的信息
    integer,allocatable             :: renumber(:)

    !! 读取 match 文件时的临时变量
    character(len=120)              :: buffer
    integer                         :: istat
    integer                         :: a,b

    !! 重新编号时的临时变量
    logical                         :: all_determined       !! 是否所有原子都已经对准
    !! atoms1相关的变量
    integer                         :: max_connect          !! 最大的链接度
    integer                         :: temp_connect         !! 保存当前原子链接度的临时变量
    integer                         :: max_connect_atom     !! 有着最大链接度的原子
    integer,allocatable             :: max_connect_atoms(:) !! 与max_connect_atom相连的原子
    !! atoms2相关的变量
    integer,allocatable             :: candidate_atoms(:)   !! atoms2 中可能与max_connect_atom对应的原子
    !! 计算二面角相关的变量
    integer                         :: a1,a2,a3     !! 计算 a1-a2-a3-max_connect_atom 的二面角
    real                            :: d_ref        !! a1-a2-a3-max_connect_atom 的二面角作为参考值
    integer                         :: b1,b2,b3     !! 计算 b1-b2-b3-candidate_atoms(i) 的二面角
    real,allocatable                :: d(:)         !! b1-b2-b3-candidate_atoms(i) 的二面角的值
    real                            :: min_def      !! 二面角的值的最小偏差
    integer                         :: f            !! 最终与max_connect_atom对应的原子

    !! 文件io单元临时变量
    integer                         :: io_in
    integer                         :: io_out

    !! 计数器
    integer                         :: i, j, k


    !! 保存命令行参数
    narg = iargc()
    if (allocated(arg)) deallocate(arg)
    allocate(arg(narg+1))
    do i=1,narg+1
        call getarg(i-1,arg(i))
    end do
    write(*,"(A)") "command line"
    write(*,"(1X,*(A,2X))") (trim(arg(i)),i=1,narg+1)
    write(*,"(A)") "-------------------------------------------------"

    !! 读取分子结构信息
    call parse_gjf(arg(2),atoms1)
    call parse_gjf(arg(3),atoms2)

    !! 输出部分分子结构信息
    write(*,"(A)") trim(arg(2))
    write(*,"(1X,I0,2X,I0)") atoms1%charge, atoms1%spin_mult
    do i=1,MIN(10,atoms1%natoms)
        write(*,"(A2,2X,3(F12.7,2X))") elemsyl(atoms1%symbols(i)), atoms1%coordinate(:,i)
    end do
    do i=1,MIN(10,atoms1%natoms)
        write(temp_string,"(A,I0)") " ",i
        conn = trim(temp_string)
        do j=i+1,atoms1%natoms
            if (atoms1%connect(i,j)>0.1) then
                write(temp_string,"(A,I0,A,F0.1)") " ",j," ",atoms1%connect(i,j)
                conn = trim(conn)//trim(temp_string)
            end if
        end do
        write(*,"(A)") trim(conn)
    end do
    write(*,"(A)") "-------------------------------------------------"
    write(*,"(A)") trim(arg(3))
    write(*,"(1X,I0,2X,I0)") atoms2%charge, atoms2%spin_mult
    do i=1,MIN(10,atoms2%natoms)
        write(*,"(A2,2X,3(F12.7,2X))") elemsyl(atoms2%symbols(i)), atoms2%coordinate(:,i)
    end do
    do i=1,MIN(10,atoms2%natoms)
        write(temp_string,"(A,I0)") " ",i
        conn = trim(temp_string)
        do j=i+1,atoms2%natoms
            if (atoms2%connect(i,j)>0.0) then
                write(temp_string,"(A,I0,A,F0.1)") " ",j," ",atoms2%connect(i,j)
                conn = trim(conn)//trim(temp_string)
            end if
        end do
        write(*,"(A)") trim(conn)
    end do
    write(*,"(A)") "-------------------------------------------------"

    !! 读取renumber_input信息
    if (allocated(renumber)) deallocate(renumber)
    allocate(renumber(atoms1%natoms))
    renumber = 0
    call open_old_file(io_in,arg(4))
    call open_new_file(io_out,"match_out")
    do
        read(io_in,"(A)",iostat=istat) buffer
        if (istat/=0) exit
        if (len_trim(buffer)==0) cycle
        read(buffer,*) a,b
        renumber(a) = b
        write(*,"(I0,2X,I0)") a,renumber(a)
        write(io_out,"(I0,2X,I0)") a,renumber(a)
    end do
    close(io_in)
    write(*,"(A)") "-------------------------------------------------"


    !! 依据链接度信息重新编号
    do

        !! 若已经完全确定，返回
        all_determined = .True.
        do i=1,size(renumber)
            if (renumber(i)==0) then
                all_determined = .False.
                exit
            end if
        end do
        if (all_determined) then
            write(*,"(A)") "match end"
            exit
        else
            write(*,"(A)") "continue matching"
        end if

        !! 找到与已经确认的原子链接最多的原子 max_connect_atom
        max_connect = -1
        do i=1,size(renumber)
            if (renumber(i)/=0) cycle !! 删除已经确认的
            temp_connect = 0    !! 初始化
            do j=1,size(renumber)
                if (i==j) cycle
                if (renumber(j)==0) cycle
                if (atoms1%connect(i,j)>0.1) temp_connect = temp_connect + 1
            end do
            if (temp_connect>max_connect) then
                max_connect      = temp_connect
                max_connect_atom = i
            end if
        end do

        write(*,"(A,I0)")       "target  atom in atoms1 : ", max_connect_atom

        if (allocated(max_connect_atoms)) deallocate(max_connect_atoms)
        do i=1,size(renumber)
            if (max_connect_atom==i)    cycle
            if (renumber(i)==0)         cycle
            if (atoms1%connect(i,max_connect_atom)>0.1) then
                if (allocated(max_connect_atoms)) then
                    max_connect_atoms = [max_connect_atoms,i]
                else
                    allocate(max_connect_atoms(1))
                    max_connect_atoms(1)= i
                end if
            end if
        end do
        if (allocated(max_connect_atoms)) then
            write(*,"(A,*(I0,1X))") "support atom in atoms1 : ", (max_connect_atoms(i),i=1,size(max_connect_atoms))
        end if

        !! 在 atoms2 中找到相应的原子们
        if (allocated(candidate_atoms)) deallocate(candidate_atoms)
        loop_i : do i=1,size(renumber)

            !! 需要与 max_connect_atom 有相同的元素符号
            if (atoms2%symbols(i)/=atoms1%symbols(max_connect_atom)) cycle loop_i

            !! 需要有相同的链接度
            if (atoms2%link(i)/=atoms1%link(max_connect_atom))       cycle loop_i

            !! 如果已经被对应了，删除
            loop_j : do j=1,size(renumber)
                if (renumber(j)==i) cycle loop_i
            end do loop_j

            !! 需要与对应后的 max_connect_atoms 有链接
            if (allocated(max_connect_atoms)) then
                do k=1,size(max_connect_atoms)
                    if (atoms2%connect(i,renumber(max_connect_atoms(k)))<0.1) cycle loop_i
                end do
            end if

            !! 接收 i
            if (allocated(candidate_atoms)) then
                candidate_atoms = [candidate_atoms,i]
            else
                allocate(candidate_atoms(1))
                candidate_atoms(1) = i
            end if
        end do loop_i

        !! 判断 candidate_atoms
        if (.not.allocated(candidate_atoms)) then
            !! candidate_atoms 为空，需要手动选择
            write(*,"(A,A,A,I0,A,A)") "Which atom in ",trim(arg(3)), " coorsponds to atom ",max_connect_atom," in ",trim(arg(2))
            read(*,*) renumber(max_connect_atom)
            write(io_out,"(I0,2X,I0)") max_connect_atom, renumber(max_connect_atom)


        else if (size(candidate_atoms)==1) then
            renumber(max_connect_atom) = candidate_atoms(1)
            write(*,"(A,I0)") "match atom in atoms2 : ", candidate_atoms(1)
            write(io_out,"(I0,2X,I0)") max_connect_atom, renumber(max_connect_atom)

        else
            !! 用二面角判断

            if (allocated(d)) deallocate(d)
            allocate(d(size(candidate_atoms)))

            !! 找到一个4个原子的序列
            loopa1 : do a1=1,size(renumber) !! 找一个与max_connect_atom相连的

                if (a1==max_connect_atom)                       cycle loopa1
                if (renumber(a1)==0)                            cycle loopa1
                if (atoms1%connect(max_connect_atom,a1)<0.1)    cycle loopa1

                loopa2 : do a2=1,size(renumber) !!找一个与a1相连的
                    if (a2==max_connect_atom)       cycle loopa2
                    if (a2==a1)                     cycle loopa2
                    if (renumber(a2)==0)            cycle loopa2
                    if (atoms1%connect(a2,a1)<0.1)  cycle loopa2

                    loopa3 : do a3=1,size(renumber) !! 找一个与a2相连的
                        if (a3==max_connect_atom)       cycle loopa3
                        if (a3==a1)                     cycle loopa3
                        if (a3==a2)                     cycle loopa3
                        if (renumber(a3)==0)            cycle loopa3
                        if (atoms1%connect(a2,a3)<0.1)  cycle loopa3
                        
                        !! 至此a1,a2,a3满足全部要求
                        write(*,"(3(A,I0,2X))") "a1=",a1,"a2=",a2,"a3=",a3
                        
                        d_ref = dihedral( atoms1%coordinate(:,a3), atoms1%coordinate(:,a2), &
                                        & atoms1%coordinate(:,a1), atoms1%coordinate(:,max_connect_atom) )

                        b1 = renumber(a1)
                        b2 = renumber(a2)
                        b3 = renumber(a3)
                        do i=1,size(candidate_atoms)
                            d(i) = dihedral( atoms2%coordinate(:,b3), atoms2%coordinate(:,b2), &
                                            & atoms2%coordinate(:,b1), atoms2%coordinate(:,candidate_atoms(i)) )
                            
                        end do

                        min_def = 100000.0
                        do i=1,size(candidate_atoms)
                            if ( mod_dihe(d(i)-d_ref) < min_def) then
                                f = candidate_atoms(i)
                                min_def = mod_dihe(d(i)-d_ref)
                            end if
                        end do
                        renumber(max_connect_atom) = f

                        write(*,"(A,I0)") "match atom in atoms2 by dihedral : ", f
                        write(io_out,"(I0,2X,I0)") max_connect_atom, renumber(max_connect_atom)

                        goto 1000

                    end do loopa3
                end do loopa2
            end do loopa1

            !! 说明不能找到a1-a2-a3
            write(*,"(A,A,A,I0,A,A)") "Which atom in ",trim(arg(3)), " coorsponds to atom ",max_connect_atom," in ",trim(arg(2))
            write(*,"(A,*(I0,1X))") "suggested value : ", (candidate_atoms(i),i=1,size(candidate_atoms))
            read(*,*) renumber(max_connect_atom)
            write(io_out,"(I0,2X,I0)") max_connect_atom, renumber(max_connect_atom)

        end if

1000    continue    !! 已经根据二面角选定了对应原子
    end do

    close(io_out)

contains

!*******************************************************************************
!> 打开一个已经存在的文件
subroutine open_old_file(io,filename)
    implicit none
    integer,intent(out)         :: io 
    character(len=*),intent(in) :: filename 
    integer                     :: istat
    character(len=80)           :: msg
    open(newunit=io,file=trim(filename),status='old',iostat=istat,iomsg=msg)
    if (istat/=0) then
        write(*,"(A,I0)") "open_old_file "//trim(filename)//" failed: istat = ", istat
        write(*,"(A)")    "error message = ", trim(msg)
        write(*,"(A)")    "shutting down..."
        stop
    end if 
    return
end subroutine open_old_file
!*******************************************************************************


!*******************************************************************************
!> 读取gjf文件的信息
subroutine parse_gjf(filename,a)
    implicit none
    character(len=*),intent(in) :: filename
    class(atoms),intent(out)    :: a
    !---------------------------------------------------------------------------
    integer                     :: io, istat
    integer                     :: nblank
    character(len=256)          :: buffer
    integer                     :: i,j
    character(len=3)            :: symbol
    integer                     :: junk2
    integer                     :: index_temp
    integer                     :: number
    real                        :: link
    call open_old_file(io,filename)
    nblank = 0
    do
        read(io,"(A)",iostat=istat) buffer
        if (istat/=0) then
            write(*,*) "parse gjf "//trim(filename)//"failed, wrong blank line."
            stop
        end if
        if (len_trim(buffer)==0) nblank = nblank + 1
        if (nblank==2) exit
    end do
    !---------------------------------------------------------------------------
    call get_gjf_natoms(filename,a%natoms)
    if (allocated(a%symbols))    deallocate(a%symbols)
    if (allocated(a%coordinate)) deallocate(a%coordinate)
    if (allocated(a%connect))    deallocate(a%connect)
    if (allocated(a%link))       deallocate(a%link)
    allocate(a%symbols(a%natoms))
    allocate(a%coordinate(3,a%natoms))
    allocate(a%connect(a%natoms,a%natoms))
    allocate(a%link(a%natoms))
    a%connect = 0.0
    a%link = 0
    !---------------------------------------------------------------------------
    !! 读取元素符号和坐标的信息
    read(io,*) a%charge, a%spin_mult
    do i=1,a%natoms
        read(io,"(A)",iostat=istat) buffer
        if ( index(buffer,"(")/=0 .and. index(buffer,")")/=0 ) then
            !! treat as PDB
            buffer = adjustl(buffer)
            index_temp = index(buffer,"(")
            symbol = buffer(1:index_temp-1)
            a%symbols(i) = elemord(symbol)

            index_temp = index(buffer,")")
            buffer = buffer(index_temp+1:)
            buffer = adjustl(buffer)
            read(buffer,*) junk2, (a%coordinate(j,i), j=1,3)
        else
            read(buffer,*) symbol, (a%coordinate(j,i), j=1,3)
            a%symbols(i) = elemord(symbol)
        end if
    end do
    !---------------------------------------------------------------------------
    !! 读取连接度信息
    read(io,*)
    do i=1,a%natoms
        read(io,"(A)",iostat=istat) buffer
        buffer = adjustl(buffer)
        index_temp = index(buffer," ")
        buffer = buffer(index_temp+1:)
        buffer = adjustl(buffer)
        do
            if (len_trim(buffer)==0) exit

            index_temp = index(buffer," ")
            read(buffer(1:index_temp-1),*) number
            buffer = buffer(index_temp+1:)
            buffer = adjustl(buffer)

            index_temp = index(buffer," ")
            read(buffer(1:index_temp-1),*) link
            buffer = buffer(index_temp+1:)
            buffer = adjustl(buffer)

            a%connect(i,number) = link
            a%connect(number,i) = link

            a%link(i)       = a%link(i) + 1
            a%link(number)  = a%link(number) + 1
        end do
    end do
    close(io)
    !---------------------------------------------------------------------------
    return
end subroutine
!*******************************************************************************

!*******************************************************************************
!>  得到gjf文件中的原子数
subroutine get_gjf_natoms(filename,natoms)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(out)         :: natoms
    !---------------------------------------------------------------------------
    integer                     :: io, istat
    integer                     :: nblank
    character(len=120)          :: buffer
    !---------------------------------------------------------------------------
    call open_old_file(io,filename)
    nblank = 0
    do
        read(io,"(A)",iostat=istat) buffer
        if (istat/=0) then
            write(*,*) "parse gjf "//trim(filename)//"failed, wrong blank lines."
            stop
        end if
        if (len_trim(buffer)==0) nblank = nblank + 1
        if (nblank==2) exit
    end do
    read(io,*) !! charge and spin line
    !---------------------------------------------------------------------------
    natoms = 0
    do
        read(io,"(A)",iostat=istat) buffer
        if (istat/=0) then
            write(*,*) "parse gjf "//trim(filename)//"failed, wrong blank lines."
            stop
        end if
        if (len_trim(buffer)==0) exit
        natoms = natoms + 1
    end do
    close(io)
    return
end subroutine get_gjf_natoms
!*******************************************************************************


function elemsyl(n)
    implicit none
    character(len=2)   :: elemsyl
    integer,intent(in) :: n
    character(len=2)   :: elem(109) =                             &
         [    'H ','He','Li','Be','B ',  'C ','N ','O ','F ','Ne',&
         &    'Na','Mg','Al','Si','P ',  'S ','Cl','Ar','K ','Ca',&
         &    'Sc','Ti','V ','Cr','Mn',  'Fe','Co','Ni','Cu','Zn',&
         &    'Ga','Ge','As','Se','Br',  'Kr','Rb','Sr','Y ','Zr',&
         &    'Nb','Mo','Tc','Ru','Rh',  'Pd','Ag','Cd','In','Sn',&
         &    'Sb','Te','I ','Xe','Cs',  'Ba','La','Ce','Pr','Nd',&
         &    'Pm','Sm','Eu','Gd','Tb',  'Dy','Ho','Er','Tm','Yb',&
         &    'Lu','Hf','Ta','W ','Re',  'Os','Ir','Pt','Au','Hg',&
         &    'Tl','Pb','Bi','Po','At',  'Rn','Fr','Ra','Ac','Th',&
         &    'Pa','U ','Np','Pu','Am',  'Cm','Bk','Cf','Es','Fm',&
         &    'Md','No','Lr','Rf','Db',  'Sg','Bh','Hs','Mt'      ]
    if (n>0.and.n<110) then
        elemsyl=elem(n)
        return
    else
        write(*,*) 'Wrong nuclear charge : ',n
        stop
    end if
end function elemsyl

function elemord(symbol)
    implicit none
    integer                     :: elemord
    character(len=*),intent(in) :: symbol
    integer                     :: i
    character(len=2)            :: symbols1(109) =                &
         [    'H ','He','Li','Be','B ',  'C ','N ','O ','F ','Ne',&
         &    'Na','Mg','Al','Si','P ',  'S ','Cl','Ar','K ','Ca',&
         &    'Sc','Ti','V ','Cr','Mn',  'Fe','Co','Ni','Cu','Zn',&
         &    'Ga','Ge','As','Se','Br',  'Kr','Rb','Sr','Y ','Zr',&
         &    'Nb','Mo','Tc','Ru','Rh',  'Pd','Ag','Cd','In','Sn',&
         &    'Sb','Te','I ','Xe','Cs',  'Ba','La','Ce','Pr','Nd',&
         &    'Pm','Sm','Eu','Gd','Tb',  'Dy','Ho','Er','Tm','Yb',&
         &    'Lu','Hf','Ta','W ','Re',  'Os','Ir','Pt','Au','Hg',&
         &    'Tl','Pb','Bi','Po','At',  'Rn','Fr','Ra','Ac','Th',&
         &    'Pa','U ','Np','Pu','Am',  'Cm','Bk','Cf','Es','Fm',&
         &    'Md','No','Lr','Rf','Db',  'Sg','Bh','Hs','Mt'      ]
    character(len=2)            :: symbols2(109) =                &
         [    'h ','he','li','be','b ',  'c ','n ','o ','f ','ne',&
         &    'na','mg','al','si','p ',  's ','cl','ar','k ','ca',&
         &    'sc','ti','v ','cr','mn',  'fe','co','ni','cu','zn',&
         &    'ga','ge','as','se','br',  'kr','rb','sr','y ','zr',&
         &    'nb','mo','tc','ru','rh',  'pd','ag','cd','in','sn',&
         &    'sb','te','i ','xe','cs',  'ba','la','ce','pr','nd',&
         &    'pm','sm','eu','gd','tb',  'dy','ho','er','tm','yb',&
         &    'lu','hf','ta','w ','re',  'os','ir','pt','au','hg',&
         &    'tl','pb','bi','po','at',  'rn','fr','ra','ac','th',&
         &    'pa','u ','np','pu','am',  'cm','bk','cf','es','fm',&
         &    'md','no','lr','rf','db',  'sg','bh','hs','mt'      ]
    elemord = 0
    do i=1,109
        if ( trim(symbol)==trim(symbols1(i)) .or. trim(symbol)==trim(symbols2(i)) ) then
            elemord = i
            exit
        end if
    end do
    if (elemord==0) then
        write(*,*) "get atomic order failed for ",symbol
    end if
end function elemord

subroutine open_new_file(io,filename)
    implicit none
    integer,intent(out)         :: io 
    character(len=*),intent(in) :: filename 
    integer                     :: istat
    character(len=80)           :: msg
    open(newunit=io,file=trim(filename),status='replace',iostat=istat,iomsg=msg)
    if (istat/=0) then
        write(*,"(A,I0)") "open_new_file "//trim(filename)//" failed: istat = ", istat
        write(*,"(A)")    "error message = "//trim(msg)
        write(*,"(A)")    "shutting down..."
        stop
    end if
    return
end subroutine open_new_file

function dihedral(a,b,c,d)
    implicit none
    real            :: dihedral
    real,intent(in) :: a(3),b(3),c(3),d(3)
    !---------------------------------------------------------------------------
    real            :: ab(3),bc(3),cd(3)
    real            :: n1(3),n2(3)
    real            :: mod_n1,mod_n2,sgn
    real,parameter  :: PI =3.1415926535898
    !---------------------------------------------------------------------------
    ab = b-a
    bc = c-b
    cd = d-c
    !---------------------------------------------------------------------------
    n1 = CROSS_PRODUCT(ab,bc)
    n2 = CROSS_PRODUCT(bc,cd)
    !---------------------------------------------------------------------------
    mod_n1 = sqrt(n1(1)**2 + n1(2)**2 + n1(3)**2)
    mod_n2 = sqrt(n2(1)**2 + n2(2)**2 + n2(3)**2)
    !---------------------------------------------------------------------------
    dihedral = DOT_PRODUCT(n1,n2)/(mod_n1*mod_n2)
    dihedral = MIN(MAX(dihedral,-1.0),1.0) ! deal with numerical error
    dihedral = 180.0*ACOS(dihedral)/PI
    sgn = DOT_PRODUCT(ab,n2)
    if(sgn < 0.0d0) dihedral = -dihedral
    !---------------------------------------------------------------------------
    return
end function
!*******************************************************************************

function mod_dihe(d)
    !! 当二面角差接近360度，修改距离360被数最近的差
    implicit none
    real :: d
    real :: mod_dihe
    mod_dihe = MIN(ABS(d),ABS(d-360.0),ABS(d-720.0),ABS(d+360.0),ABS(d+720.0))
end function

function CROSS_PRODUCT(v1,v2)
    implicit none
    real            :: CROSS_PRODUCT(3)
    real,intent(in) :: v1(3), v2(3)
    CROSS_PRODUCT(1) =   v1(2)*v2(3)-v1(3)*v2(2)
    CROSS_PRODUCT(2) = -(v1(1)*v2(3)-v1(3)*v2(1))
    CROSS_PRODUCT(3) =   v1(1)*v2(2)-v1(2)*v2(1)
end function

end program match_gjf