
module points

    implicit none

    type :: point
        integer :: indice
        real(kind=8) :: x
        real(kind=8) :: y
    end type point

contains

    subroutine construction_points(p, xx, yy, indice)
        type(point), intent(out) :: p
        real(kind=8), intent(in) :: xx, yy
        integer, intent(in) :: indice
        p%indice = indice
        p%x = xx
        p%y = yy
    end subroutine construction_points

    function comparaison_points(pt1, pt2) result(res)
        type(point), intent(in) :: pt1, pt2
        logical :: res
        res = .false.
        if (pt1%indice == pt2%indice .and. pt1%x == pt2%x .and. pt1%y == pt2%y) then
            res = .true.
        end if
    end function comparaison_points


    subroutine nuage(xmax,xmin,ymax,ymin,n,lespoints)
        real(8),intent(in)::xmax,ymax,ymin,xmin
        integer,intent(in)::n
        type(point),dimension(:),intent(out)::lespoints
        integer ::  i,j
        real(8) :: dx,dy,x1,x2,y1,y2
        real(8), allocatable :: x(:), y(:),x_sub(:),y_sub(:)
        real(8) :: rand_x, rand_y

        x1=xmin
        x2=xmax
        y1=ymin
        y2=ymax
        call random_seed()
        allocate(x_sub(n+1),y_sub(n+1))
        x_sub(1)=x1
        y_sub(1)=y1
        dy=(y2-y1)/n
        dx=(x2-x1)/n
        do i=1,n
            x_sub(i+1)=x1+i*dx
            y_sub(i+1)=y1+i*dy
        end do
        allocate(x(n*n), y(n*n))
        do i=1,n
            do j=1,n
                call random_number(rand_x)
                x((i-1)*n+j) = x_sub(i) + rand_x * dx
                call random_number(rand_y)
                y((i-1)*n+j) = y_sub(j) + rand_y * dy
                if(j>1) then
                    do while (abs(y(j) - y(j-1)) < 1.0e-4)
                        call random_number(rand_y)
                        y(i) = y_sub(i) + rand_y * dy
                    end do
                end if
            end do
            if (i>1) then 
                do while (abs(x(i) - x(i-1)) < 1.0e-4)
                    call random_number(rand_y)
                    y(i) = y_sub(i) + rand_y * dy
                end do
            end if
        end do 
        do i = 1, n*n
            call construction_points(lespoints(i+3), x(i), y(i),i+3)
        end do
        do i = 1, 3
            call construction_points(lespoints(i), x(i), y(i),i)
        end do
        deallocate(x, y,x_sub,y_sub)
    end subroutine nuage

end module points
