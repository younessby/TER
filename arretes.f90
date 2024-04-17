module arretes
    use points
    implicit none

    type, public :: edge
        type(point) :: p1     ! Point de départ
        type(point) :: p2     ! Point final
        integer :: indice     ! Indice de l'arête
    end type edge

contains

    subroutine construction_edge(arrete, pt1, pt2, indice)
        type(edge), intent(out) :: arrete
        type(point), intent(in) :: pt1, pt2
        integer, intent(in) :: indice
        arrete%p1 = pt1
        arrete%p2 = pt2
        arrete%indice = indice
    end subroutine construction_edge

    function comparaison_arrete(arrete1, arrete2) result(res)
        type(edge), intent(in) :: arrete1, arrete2
        logical :: res
        res = .false.
        if ((comparaison_points(arrete1%p1, arrete2%p1) .and. comparaison_points(arrete1%p2, arrete2%p2)) .or. &
            (comparaison_points(arrete1%p2, arrete2%p1) .and. comparaison_points(arrete1%p1, arrete2%p2))) then
            res = .true.
        end if
    end function comparaison_arrete

    function comparaison_edge(arrete1, arrete2) result(res)
        type(edge), intent(in) :: arrete1, arrete2
        logical :: res
        res = .false.
        if (arrete1%indice == arrete2%indice) then
            res = .true.
        end if
    end function comparaison_edge

end module arretes
