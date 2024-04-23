module triangulation

    use points
    use arretes
    use triangles
    
    implicit none
    
    contains

    
    Function SuperTriangle(vertices) result(superTri)
        type(point), intent(in), dimension(:) :: vertices ! les vertices sont des points
        type(triangle) :: superTri
        real(kind=8) :: xmin, xmax, ymin, ymax
        integer :: i
        xmin = vertices(1)%x
        xmax = vertices(1)%x
        ymin = vertices(1)%y
        ymax = vertices(1)%y
        do i = 2, size(vertices)
            xmin = min(xmin, vertices(i)%x)
            xmax = max(xmax, vertices(i)%x)
            ymin = min(ymin, vertices(i)%y)
            ymax = max(ymax, vertices(i)%y)
        end do
        superTri%p1%x = xmin - 1.0
        superTri%p1%y = ymin - 1.0
        superTri%p1%indice = 1
        superTri%p2%x = xmin + 2.0 * (xmax - xmin) + 3.0
        superTri%p2%y = ymin - 1.0
        superTri%p2%indice = 2
        superTri%p3%x = xmin - 1.0
        superTri%p3%y = ymin + 2.0 * (ymax - ymin) + 3.0
        superTri%p3%indice = 3
        superTri%indice = 0
    end Function SuperTriangle
    
    
subroutine algo_tri(pts,tris,n)
    integer,intent(in)::n
    type(point),dimension(:),intent(in)::pts
    type(triangle),dimension(:),intent(out)::tris
    type(triangle),allocatable,dimension(:)::bad_trianagle,tri_0,tri_1,nett
    type(edge),allocatable,dimension(:)::all_edge,polygone,tab_tri
    type(triangle)::super,new_tri,new_tri1
    integer::i,j,nmtri,nmedge,nmpol,nmbad,jj,nmtr0,nett0,nm
    type(point):: ppoint
    type(triangle):: tri0 ,supertri

    nmtri=1
    tri0=tris(1)
    supertri=tris(1)
    do i=4,n
        ppoint=pts(i)
        allocate(polygone(3*n),bad_trianagle(3*n),nett(3*n))
        do j=1,3*n
            nett(i)=supertri
        end do
        nmbad=0
        nmtr0=nmtri
        nett0=0
        do j=1,nmtri
            if (CircumCircleContains(tris(j),ppoint)) then
                nmbad=nmbad+1
                bad_trianagle(nmbad)=tris(j)
                nmtr0=nmtr0-1
            else 
                nett0=nett0+1
                tris(j)%indice=nett0
                nett(nett0)=tris(j)
            end if
        end do 
        nmtri=nmtr0
        allocate(all_edge(3*n))
        nmedge=0
        do j=1,nmbad
            nmedge=nmedge+1
             call construction_edge(all_edge(nmedge),bad_trianagle(j)%p1,bad_trianagle(j)%p2,nmedge)
            nmedge=nmedge+1
            call construction_edge(all_edge(nmedge),bad_trianagle(j)%p3,bad_trianagle(j)%p2,nmedge)
            nmedge=nmedge+1
            call construction_edge(all_edge(nmedge),bad_trianagle(j)%p1,bad_trianagle(j)%p3,nmedge)
        end do 
        nmpol=0
        do j=1,nmedge
            if (count(all_edge(1:nmedge),all_edge(j),nmedge)==1) then 
                nmpol=nmpol+1
                polygone(nmpol)=all_edge(j)
            end if 
        end do 
        tris(1:nett0)=nett(1:nett0)
        do j=1,nmpol
            nett0=nett0+1
            call construction_triangle(new_tri,polygone(j)%p1,polygone(j)%p2,ppoint,nett0)
            call sens_tri(new_tri,new_tri1)
            tris(nett0)=new_tri1
        end do
        nmtri=nett0
        deallocate(polygone,all_edge,bad_trianagle,nett)
    end do
end subroutine algo_tri


function count(array, element,n) result(cnt)
    implicit none
    integer,intent(in)::n
    type(edge), intent(in) :: array(:)
    type(edge), intent(in) :: element
    integer :: i, cnt

    cnt = 0
    do i = 1, n
        if (comparaison_arrete(array(i), element)) then
            cnt = cnt + 1
        end if
    end do
end function count

subroutine trier_et_supprimer_doublons(array, n,nm) 
    implicit none
    integer, intent(in) :: n
    integer,intent(out)::nm
    type(edge),dimension(:), intent(inout) :: array
    type(edge),dimension(:),allocatable :: tab
    integer :: i, j, k,res
    allocate(tab(n))
    nm=1
    tab(1)=array(1)
    do i = 2, n
        res=0
        do j=1,nm
            if (comparaison_arrete(array(i), tab(j))) then
                res=1
            end if
        end do 
        if (res==0) then 
            nm=nm+1
            tab(nm) = array(i)
        end if 
    end do
    array(1:nm)=tab(1:nm)
    deallocate(tab)
end subroutine trier_et_supprimer_doublons

subroutine connectivite(lespoints,lestriangles,lesaretes,aretes_maille,maille_aretes,maille_point,nb_point,nb_tri,nb_arete)
    type(point),dimension(:),intent(in)::lespoints
    type(triangle),dimension(:),intent(in)::lestriangles
    type(edge),dimension(:), intent(out)::lesaretes
    integer, dimension(:,:), intent(out) :: aretes_maille,maille_aretes
    integer,dimension(:,:),intent(inout)::maille_point
    integer,intent(in)::nb_tri,nb_point
    integer,intent(out)::nb_arete
    type(edge),dimension(3)::edges1,edges2

    integer::i,nb_edge,pp1,pp2,pp3,nbtri_2,nbtri_1,j
    integer, dimension(:,:), allocatable :: tab
    allocate(tab(1:nb_point,nb_point))
    tab=0
    nb_edge=0
    do i=1,nb_tri
        pp1=lestriangles(i)%p1%indice
        pp2=lestriangles(i)%p2%indice
        pp3=lestriangles(i)%p3%indice
        maille_point(pp1,1)=maille_point(pp1,1)+1
        maille_point(pp1,maille_point(pp1,1)+1)=i
        maille_point(pp2,1)=maille_point(pp2,1)+1
        maille_point(pp2,maille_point(pp2,1)+1)=i
        maille_point(pp3,1)=maille_point(pp3,1)+1
        maille_point(pp3,maille_point(pp3,1)+1)=i
    end do 


    do i=1,nb_tri
        pp1=lestriangles(i)%p1%indice
        pp2=lestriangles(i)%p2%indice
        pp3=lestriangles(i)%p3%indice
        if(tab(pp1,pp2)==0 .and. tab(pp2,pp1)==0) then
            nb_edge=nb_edge+1
            tab(pp1,pp2)=nb_edge
            tab(pp2,pp1)=nb_edge
            call construction_edge(lesaretes(nb_edge),lestriangles(i)%p1,lestriangles(i)%p2,nb_edge)
        end if
        if(tab(pp3,pp2)==0 .and. tab(pp2,pp3)==0) then
            nb_edge=nb_edge+1
            tab(pp3,pp2)=nb_edge
            tab(pp2,pp3)=nb_edge
            call construction_edge(lesaretes(nb_edge),lestriangles(i)%p3,lestriangles(i)%p2,nb_edge)
        end if
        if(tab(pp3,pp1)==0 .and. tab(pp1,pp3)==0) then
            nb_edge=nb_edge+1
            tab(pp3,pp1)=nb_edge
            tab(pp1,pp3)=nb_edge
            call construction_edge(lesaretes(nb_edge),lestriangles(i)%p1,lestriangles(i)%p3,nb_edge)
        end if
    end do
    nb_arete=nb_edge
    !allocate(aretes_maille(1:nb_tri,3),maille_aretes(1:nb_edge,2))
        aretes_maille=0
        maille_aretes=0
        do i=1,nb_tri
            pp1=lestriangles(i)%p1%indice
            pp2=lestriangles(i)%p2%indice
            pp3=lestriangles(i)%p3%indice

            aretes_maille(i,1)=tab(pp1,pp2)
            if(maille_aretes(aretes_maille(i,1),1)==0) then 
                maille_aretes(aretes_maille(i,1),1)=i
            else
                maille_aretes(aretes_maille(i,1),2)=i
            end if 

            aretes_maille(i,2)=tab(pp2,pp3)
            if(maille_aretes(aretes_maille(i,2),1)==0) then 
                maille_aretes(aretes_maille(i,2),1)=i
            else
                maille_aretes(aretes_maille(i,2),2)=i
            end if 

            aretes_maille(i,3)=tab(pp3,pp1)
            if(maille_aretes(aretes_maille(i,3),1)==0) then 
                maille_aretes(aretes_maille(i,3),1)=i
            else
                maille_aretes(aretes_maille(i,3),2)=i
            end if 
        end do
        do i=1,nb_arete
            nbtri_1=maille_aretes(i,1)
            nbtri_2=maille_aretes(i,2)
            if(arret_sens(lestriangles(nbtri_1), lesaretes(i)).eqv. .false.) then 
                maille_aretes(i,1)=nbtri_2
                maille_aretes(i,2)=nbtri_1
            end if
        end do 



        deallocate(tab)

    end subroutine connectivite

    subroutine construction_points_center(domingo, p_center)
        type(triangle), intent(in)     :: domingo
        type(point), intent(out)       :: p_center
        real(8)                        :: x1, x2, x3, y1, y2, y3

        x1 = domingo%p1%x
        x2 = domingo%p2%x
        x3 = domingo%p3%x
        y1 = domingo%p1%y
        y2 = domingo%p2%y
        y3 = domingo%p3%y

        p_center%x=((x1**2+y1**2)*(y2-y3)+(x2**2+y2**2)*(y3-y1)+(x3**2+y3**2)*(y1-y2))/(2.0*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)))
        p_center%y=((x1**2+y1**2)*(x3-x2)+(x2**2+y2**2)*(x1-x3)+(x3**2+y3**2)*(x2-x1))/(2.0*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)))
        p_center%indice=domingo%indice

    end subroutine










            










end module triangulation