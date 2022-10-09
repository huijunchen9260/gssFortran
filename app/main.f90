program main

    use iso_Fortran_env, only: rk => real64, ik => int32
    use goldenSectionSearch
    implicit none
    integer(ik), parameter :: bknum = 100
    integer(ik), parameter :: knum = 100
    integer(ik), parameter :: enum = 100
    real(rk) :: fout, xout, t0, t1



    type configurations
        real(rk) :: addition
        real(rk) :: zval                    ! TFP shock
        real(rk) :: wval                    ! wage
        real(rk) :: pval                    ! Marginal utility
        real(rk) :: pfval                   ! future Marginal utility
        real(rk) :: qbval                   ! bond price
        real(rk) :: Qbuy                    ! capital purchasing price
        real(rk) :: qsell                   ! capital selling price
        real(rk) :: xuval                   ! cash on hand for upward-adj firm
        real(rk) :: xdval                   ! cash on hand for downward-adj firm
        ! real(rk) :: ewbk(bknum, knum)
        real(rk), allocatable :: ewk(:), ewbk(:, :), evbk(:, :), evk(:)
        ! real(rk) :: evbk(bknum, knum)
        ! real(rk) :: evk(knum)
    end type

    type(configurations) :: conf

    allocate(conf%ewk(knum), conf%ewbk(bknum, knum))

    conf%addition = 2.0_rk
    conf%qsell = 0.918790087890625_rk
    conf%wval = 1.087975585937500_rk

    deallocate(conf%ewk)

    call cpu_time(t0)
    call goldenSec(fout, xout, LB = -2.0_rk, UB = 7.0_rk, conf = conf)
    call cpu_time(t1)
    write(*, *) fout, xout
    write(*, *) "subroutine gss took ", t1 - t0, " seconds"

    call cpu_time(t0)
    call gss(fout, xout, f, LB = -2.0_rk, UB = 7.0_rk, show = .false., func_data = conf)
    call cpu_time(t1)
    write(*, *) fout, xout
    write(*, *) "gss with class took ", t1 - t0, " seconds"


contains

    subroutine goldenSec(fout, xout, LB, UB, conf)
        real(rk), parameter :: goldenRatio = (3.0_rk - dsqrt(5.0_rk)) / 2.0_rk
        real(rk), intent(out) :: fout, xout
        real(rk), intent(in) :: LB, UB
        type(configurations), intent(in), optional :: conf
        integer(ik) :: iter, gssMaxIter
        logical :: isFindMax, isShow
        character(len=128) :: tmp
        real(rk) :: a, b, c, d, z, fval, fc, fd, gssTol, gssDist

        iter = 0
        gssMaxIter = 500
        gssTol = 1D-9
        gssDist = 2.0_rk*gssTol
        isFindMax = .true.
        isShow = .false.

        a = LB
        b = UB
        c = a + goldenRatio*(b-a);
        d = a + (1-goldenRatio)*(b-a);

        fval = -c**2 + conf%addition; fc = -fval
        fval = -d**2 + conf%addition; fd = -fval

        do while (gssDist > gssTol .and. iter <= gssMaxIter)

            iter = iter + 1

            if (fc >= fd) then
                z = c + (1.0_rk-goldenRatio)*(b-c)
                ! case 1 [a c d b] <--- [c d z b]
                a = c
                c = d
                fc = fd
                d = z
                ! fval = Obj(d, func_data); fd = merge(-fval, fval, isFindMax);
                fval = -d**2 + conf%addition; fd = -fval
            else
                z = a + goldenRatio*(d-a)
                ! case 2 [a c d b] <--- [a z c d]
                b = d
                d = c
                fd = fc
                c = z
                fval = -c**2 + conf%addition; fc = -fval
            endif

            gssDist = b - a

            if (iter == 1 .and. isShow) then
                write(*, '(a4, 5(a20))') 'iter', '||b-a||', 'a', 'c', 'd', 'b'
                write(tmp, '(a4, 5(a20))') 'iter', '||b-a||', 'a', 'c', 'd', 'b'
                write(*, '(a)') repeat('=', len_trim(tmp))
            endif
            if (isShow) write (*, '(i4, 5(ES20.6))') iter, gssDist, a, c, d, b

        enddo

        fout = fval
        xout = z


    end subroutine goldenSec

    function f(k, func_data)
        real(rk), intent(in) :: k
        class(*), intent(in), target, optional :: func_data
        real(rk) :: f
        type(configurations) :: conf

        select type(func_data)
            type is (configurations)
            conf = func_data
        end select

        f = -k**2 + conf%addition
    end function f


end program main
