module triangles
    use points
    use arretes
    implicit none 

    type, public :: triangle
        type(point) :: p1
        type(point) :: p2    ! Les indices des sommets 
        type(point) :: p3
        integer :: indice  ! L'indice du triangle
    end type triangle


contains

    subroutine construction_triangle(triangle0, pt1, pt2, pt3, indice0)
        type(triangle), intent(out) :: triangle0
        type(point), intent(in) :: pt1, pt2, pt3
        integer, intent(in) :: indice0
        triangle0%p1 = pt1
        triangle0%p2 = pt2
        triangle0%p3 = pt3
        triangle0%indice = indice0
    end subroutine construction_triangle

    function CircumCircleContains(tri, pt) result(res)
        type(triangle), intent(in) :: tri
        type(point), intent(in) :: pt
        real(kind=8) :: x1, y1, x2, y2, x3, y3
        real(kind=8) :: D, Ux, Uy, radius, dis
        logical :: res

        x1 = tri%p1%x
        y1 = tri%p1%y
        x2 = tri%p2%x
        y2 = tri%p2%y
        x3 = tri%p3%x
        y3 = tri%p3%y

        D = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        if (abs(D) < 1.0E-10) then
            print*, "Les points sont colinéaires ou le dénominateur est trop proche de zéro."
            res = .false.
            return
        endif
        ! Calcul des coordonnées du centre du cercle circonscrit
        Ux = ((x1*x1 + y1*y1) * (y2 - y3) + (x2*x2 + y2*y2) * (y3 - y1) &
        + (x3*x3 + y3*y3) * (y1 - y2)) / D
        Uy = ((x1*x1 + y1*y1) * (x3 - x2) + (x2*x2 + y2*y2) * (x1 - x3)&
         + (x3*x3 + y3*y3) * (x2 - x1)) / D
        ! Calcul du rayon du cercle circonscrit
        radius = sqrt(((x2 - x1)**2 + (y2 - y1)**2) * ((x3 - x2)**2 &
        + (y3 - y2)**2) * ((x1 - x3)**2 + (y1 - y3)**2)) / abs(D)
        ! Calcul de la distance entre le centre du cercle circonscrit et le point
        dis = sqrt((Ux - pt%x)**2 + (Uy - pt%y)**2)
        if (dis < radius) then
            res = .true.
        else 
            res = .false.
        end if 
    end function CircumCircleContains

    subroutine sens_tri(tri1,tri2)
        type(triangle),intent(in)::tri1
        type(triangle),intent(out)::tri2
        real(8),dimension(3):: v1,v2,v3,v_ref
        integer::i
        real(8)::a
        
        v1(1)=tri1%p1%x-tri1%p2%x
        v2(1)=tri1%p1%x-tri1%p3%x
        v1(2)=tri1%p1%y-tri1%p2%y
        v2(2)=tri1%p1%y-tri1%p3%y
        v1(3)=0.0
        v2(3)=0.0

        v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
        v3(2)=-(v1(1)*v2(3)-v1(3)*v2(1))
        v3(3)=v1(1)*v2(2)-v1(2)*v1(1)

        v_ref(1)=0.0
        v_ref(2)=0.0
        v_ref(3)=0.0

        a=v3(1)*v_ref(1)+v3(2)*v_ref(2)+v3(3)*v_ref(3)

        if (a>0) then
            tri2=tri1
        else 
            tri2%p1=tri1%p1
            tri2%p2=tri1%p3
            tri2%p3=tri1%p2
        end if 
    
    end subroutine 
    function sens_tri2(tri1) result(logic)
        type(triangle),intent(in)::tri1
        real(8),dimension(3):: v1,v2,v3,v_ref
        logical::logic
        integer::i
        real(8)::a
        
        v1(1)=tri1%p1%x-tri1%p2%x
        v2(1)=tri1%p1%x-tri1%p3%x
        v1(2)=tri1%p1%y-tri1%p2%y
        v2(2)=tri1%p1%y-tri1%p3%y
        v1(3)=0.0
        v2(3)=0.0

        v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
        v3(2)=-(v1(1)*v2(3)-v1(3)*v2(1))
        v3(3)=v1(1)*v2(2)-v1(2)*v1(1)

        v_ref(1)=0.0
        v_ref(2)=0.0
        v_ref(3)=0.0

        a=v3(1)*v_ref(1)+v3(2)*v_ref(2)+v3(3)*v_ref(3)

        if (a>0) then
            logic=.true.
        else 
            logic=.false.
        end if 
    
    end function sens_tri2 

    function arret_sens(tri, arretee) result(logic)
        type(triangle),intent(in)::tri
        type(edge),intent(in)::arretee
        logical::logic
        type(point):: p1,p2,p3,pp1,pp2

        p1=tri%p1
        p2=tri%p2
        p3=tri%p3
        pp1=arretee%p1
        pp2=arretee%p2
        logic=.false.
        if (p1%indice==pp1%indice .and. p2%indice==pp2%indice) then 
            logic=.true.
        end if 
        if (p2%indice==pp1%indice .and. p3%indice==pp2%indice) then 
            logic=.true.
        end if 
        if (p3%indice==pp1%indice .and. p3%indice==pp2%indice) then 
            logic=.true.
        end if 

    end function

    subroutine construire_new(tri,new_tri,points)
        type(point),intent(in)::points
        type(triangle),intent(in):: tri
        type(triangle),intent(out)::new_tri


        if(tri%p1%indice==points%indice) then
            new_tri%indice=tri%indice
            new_tri%p1=tri%p1
            new_tri%p2=tri%p2
            new_tri%p3=tri%p3
        end if  

        if(tri%p2%indice==points%indice) then
            new_tri%indice=tri%indice
            new_tri%p1=tri%p2
            new_tri%p2=tri%p3
            new_tri%p3=tri%p1
        end if  

        if(tri%p3%indice==points%indice) then
            new_tri%indice=tri%indice
            new_tri%p1=tri%p3
            new_tri%p2=tri%p1
            new_tri%p3=tri%p2
        end if  
    end subroutine construire_new

    subroutine trouver_ar(tri,arret,pointd,edge1,edge2) 
        type(triangle),intent(in)::tri
        type(edge),dimension(:),intent(in)::arret
        type(edge),intent(out):: edge1,edge2
        type(point),intent(in)::pointd
        type(edge)::artt,artt2
        type(triangle)::new_tri
        integer::i,j

        j=0
        do i=1,3
            artt=arret(i)
            print*,"je suis dans subroutine ",pointd%indice,artt%indice,artt%indice
            if(artt%p1%indice==pointd%indice .or. artt%p2%indice==pointd%indice) then
                j=j+1 
                if(j==1) then 
                    edge1%p1=artt%p1
                    edge1%p2=artt%p2
                    edge1%indice=artt%indice

                end if
                if(j==2) then 
                    edge2%p1=artt%p1
                    edge2%p2=artt%p2
                    edge2%indice=artt%indice
                end if 
            end if 
        end do
        if (j==0) then 
            print*,"hhhhhhhhhhhhhhhhhhhhhhh"
        end if 

        if (edge1%p1%indice/=pointd%indice) then
            if(edge1%p2%indice/=pointd%indice)then
                print*,'errreru des arrt trouver 205'
                print*,edge1%p1%indice,edge1%p2%indice
                print*,edge2%p1%indice,edge2%p2%indice
                print*,pointd%indice
            end if 
            artt=edge1
            edge1=edge2
            edge2=artt
        end if 
    end subroutine trouver_ar



    ! subroutine trouver_ar2(tri,arret,edge1,pointd,edge2) 
    !     type(triangle),intent(in)::tri
    !     type(edge),dimension(:),intent(in)::arret
    !     type(edge),intent(out):: edge2
    !     type(edge),intent(in):: edge1
    !     type(point),intent(in)::pointd
    !     type(edge)::artt
    !     integer::i,j

    !     j=0
    !     do i=1,3
    !         artt=arret(i)
    !         if(artt%p1%indice==pointd%indice .or. artt%p2%indice==pointd%indice) then
    !             if (artt%indice/=edge1%indice) then 
    !                 edge2=artt
    !             end if 
    !         end if 
    !     end do
    ! end subroutine trouver_ar2

    
end module triangles
