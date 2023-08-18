module test_lines
    use constants, only: pr
    real(pr) :: self_x(374), self_y(374)
contains

    subroutine read_selfxy
        integer :: i, funit
        open(newunit=funit, file="test/self_cross_line")
        do i=1,374
            read(funit, *) self_x(i), self_y(i)
        end do
        close(funit)
    end subroutine
end module

program test_intersect
    use constants, only: pr
    use test_lines, only: self_x, self_y, read_selfxy
    use linalg, only: intersection
    implicit none
    integer, parameter :: n=2001
    real(pr) :: l1_x(n), l2_x(n)
    real(pr) :: l1_y(n), l2_y(n)
    integer :: i

    real(pr) :: inter

    l1_x = [(real(i, pr)/100._pr, i=-1000,1000)]
    l2_x = [(real(i, pr)/100._pr, i=-1000,1000)]

    l1_y = 2 * l1_x
    l2_y = l2_x ** 2

    call intersection(l1_x, l1_y, l2_x, l2_y, inter)
    call read_selfxy
    call intersection(self_x, self_y, inter)
end program

