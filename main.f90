program main

    use iso_Fortran_env, only: rk => real64, ik => int32
    use goldenSectionSearch
    implicit none
    real(rk) :: fout, xout, t0, t1

    call cpu_time(t0)
    call goldenSec(fout, xout, LB = -2.0_rk, UB = 7.0_rk)
    call cpu_time(t1)
    write(*, *) "subroutine gss took ", t1 - t0

    call cpu_time(t0)
    call gss(fout, xout, f, LB = -2.0_rk, UB = 7.0_rk, show = .false.)
    call cpu_time(t1)
    write(*, *) "gss with class took ", t1 - t0


contains

    subroutine goldenSec(fout, xout, LB, UB)
        real(rk), parameter :: goldenRatio = (3.0_rk - dsqrt(5.0_rk)) / 2.0_rk
        real(rk), intent(out) :: fout, xout
        real(rk), intent(in) :: LB, UB
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

        fval = -c**2; fc = -fval
        fval = -d**2; fd = -fval

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
                fval = -d**2; fd = -fval
            else
                z = a + goldenRatio*(d-a)
                ! case 2 [a c d b] <--- [a z c d]
                b = d
                d = c
                fd = fc
                c = z
                fval = -c**2; fc = -fval
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

        f = -k**2
    end function f


end program main
