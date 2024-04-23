
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



    subroutine removeDuplicates(arr, unique_arr, new_size)
        integer, intent(in) :: arr(:)        ! Tableau d'entrée
        integer, allocatable, intent(out) :: unique_arr(:)  ! Tableau de sortie sans doublons
        integer, intent(out) :: new_size     ! Nouvelle taille du tableau de sortie
    
        integer :: i, j
        logical :: flag                      ! Indicateur logique pour les doublons
        integer, allocatable :: temp_arr(:)
    
        ! Allocation initiale avec la taille du tableau d'entrée
        allocate(temp_arr(size(arr)))
        new_size = 0
    
        ! Boucle pour vérifier chaque élément
        do i = 1, size(arr)
          flag = .true.
          ! Vérifier s'il est déjà dans temp_arr
          do j = 1, new_size
            if (arr(i) == temp_arr(j)) then
              flag = .false.
              exit
            endif
          enddo
          ! Si l'élément n'est pas un doublon, l'ajouter à temp_arr
          if (flag) then
            new_size = new_size + 1
            temp_arr(new_size) = arr(i)
          endif
        enddo
    
        ! Redimensionner le tableau unique_arr pour contenir uniquement les éléments uniques
        allocate(unique_arr(new_size))
        unique_arr = temp_arr(1:new_size)
    
        ! Libération de la mémoire temporaire
        deallocate(temp_arr)
      end subroutine removeDuplicates



      subroutine sortPointsCCW(center, points, sorted_points)
        type(point), intent(in) :: center           ! Centre pré-déterminé
        type(point), intent(in) :: points(:)        ! Tableau de points
        type(point), intent(out) :: sorted_points(:) ! Tableau des points triés
        integer :: i, j, n
        real(kind=8), allocatable :: angles(:)
        integer, allocatable :: indices(:)
    
        n = size(points)
        allocate(angles(n))
        allocate(indices(n))
    
        ! Calcul des angles pour chaque point par rapport au centre
        do i = 1, n
          angles(i) = atan2(points(i)%y - center%y, points(i)%x - center%x)
          indices(i) = i
        end do
    
        ! Tri des indices basés sur les angles avec un tri simple (à bulles)
        call bubbleSort(angles, indices, n)
    
        ! Réorganisation des points selon les indices triés
        do i = 1, n
          sorted_points(i) = points(indices(i))
        end do
      end subroutine sortPointsCCW
    
      ! Procédure pour trier les indices selon les angles avec tri à bulles
      subroutine bubbleSort(angles, indices, n)
        real(kind=8), intent(inout) :: angles(:)
        integer, intent(inout) :: indices(:)
        integer, intent(in) :: n
        integer :: i, j
        real(kind=8) :: temp_angle
        integer :: temp_index
    
        do i = n-1, 1, -1
          do j = 1, i
            if (angles(j) > angles(j+1)) then
              ! Échanger les angles
              temp_angle = angles(j)
              angles(j) = angles(j+1)
              angles(j+1) = temp_angle
              ! Échanger les indices
              temp_index = indices(j)
              indices(j) = indices(j+1)
              indices(j+1) = temp_index
            endif
          end do
        end do
      end subroutine bubbleSort

      function hasObtuseAngle(p1, p2, p3) result(isObtuse)
        type(point), intent(in) :: p1, p2, p3
        logical :: isObtuse
        real(kind=8) :: a, b, c
    
        ! Calcul des longueurs des côtés du triangle à partir des points
        a = distance(p2, p3)
        b = distance(p1, p3)
        c = distance(p1, p2)
    
        ! Vérification pour chaque configuration de côtés
        isObtuse = (a**2 > b**2 + c**2) .or. (b**2 > a**2 + c**2) .or. (c**2 > a**2 + b**2)
      end function hasObtuseAngle
    
      function distance(p1, p2) result(dist)
        type(point), intent(in) :: p1, p2
        real(kind=8) :: dist
    
        dist = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2)
      end function distance

      function findMidpoint(p1, p2) result(mid)
        type(point), intent(in) :: p1, p2
        type(point) :: mid
    
        mid%x = 0.5 * (p1%x + p2%x)
        mid%y = 0.5 * (p1%y + p2%y)
      end function findMidpoint
    
      function getMidpointOfObtuseSide(p1, p2, p3) result(midpoint)
        type(point), intent(in) :: p1, p2, p3
        type(point) :: midpoint
        real(kind=8) :: a, b, c
    
        a = distance(p2, p3)
        b = distance(p1, p3)
        c = distance(p1, p2)
    
        if (c**2 > a**2 + b**2) then
          midpoint = findMidpoint(p2, p3)
        elseif (b**2 > a**2 + c**2) then
          midpoint = findMidpoint(p1, p3)
        else
          midpoint = findMidpoint(p1, p2)
        endif
      end function getMidpointOfObtuseSide
end module points
