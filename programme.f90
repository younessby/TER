program programme

    use points
    use arretes
    use triangles
    use triangulation

    implicit none
    type(triangle) :: super, tria
    integer :: i, n , nb_tri,np, nb_tri1,nb_aretes,nb_point,nb_edg0,j,nbr_tri_point,nb_voro,indice_tri
    real (8):: xmin, xmax, ymin, ymax ! Bornes de l'intervalle
    real(8) :: dx, dy ! Plage de variation en x et y
    type(point),allocatable,dimension(:)::lespoints,lescentres
    type(edge),allocatable,dimension(:)::lesaretes
    type(triangle),allocatable,dimension(:)::lestriangles,triangles_tri_super,triangles_tri
    integer,allocatable, dimension(:,:):: aretes_maille,maille_aretes,voronoi
    integer,allocatable,dimension(:,:)::maille_point
    type(edge)::aretes_deb,arretes_ari
    type(point)::p1,p2,p3
    type(edge)::edge1,edge2,edge3
    type(triangle)::newtri
    type(edge),dimension(3)::edges



    
    !nombre des points pour chaque intervalle 
    np=5
    n=np*np
    nb_point=n+3
    ! Bornes de l'intervalle
    xmin = -5.0
    xmax = 5.0
    ymin = -5.0
    ymax = 5.0

    allocate(lespoints(nb_point))
    !initialiser la nuage de point
    call nuage(xmax, xmin, ymax, ymin, np, lespoints)
    print*,"c'est bon"
    super = SuperTriangle(lespoints)
    call construction_points(lespoints(1),super%p1%x,super%p1%y,1)
    call construction_points(lespoints(2),super%p2%x,super%p2%y,2)
    call construction_points(lespoints(3),super%p3%x,super%p3%y,3)
    
    allocate(lestriangles(3*nb_point))

    do i=1,3*nb_point
        lestriangles(i)=super
    end do

    call algo_tri(lespoints,lestriangles,nb_point)

    !trier les triangles 
    allocate(triangles_tri_super(3*nb_point))
    nb_tri=1
    triangles_tri_super(1)=super
    do i=1,3*nb_point
       if(lestriangles(i)%p1%indice/=1 .or. lestriangles(i)%p2%indice/=2 &
       .or. lestriangles(i)%p3%indice/=3)then
        triangles_tri_super(i+1)=lestriangles(i)
            nb_tri=nb_tri+1
       end if 
    end do 

    !call save_vtk(triangles_tri_super(1:nb_tri),lespoints,"testt.vtk",n+3,nb_tri)

    allocate(triangles_tri(3*nb_tri))
        nb_tri1=0
        do i=1,nb_tri
            if(triangles_tri_super(i)%p1%indice/=1 .and. triangles_tri_super(i)%p1%indice/=2 &
            .and. triangles_tri_super(i)%p1%indice/=3 .and. triangles_tri_super(i)%p2%indice/=1 &
            .and. triangles_tri_super(i)%p2%indice/=2 .and. triangles_tri_super(i)%p2%indice/=3 &
            .and. triangles_tri_super(i)%p3%indice/=1 .and. triangles_tri_super(i)%p3%indice/=2 &
            .and. triangles_tri_super(i)%p3%indice/=3)then
                nb_tri1=nb_tri1+1
                triangles_tri(nb_tri1)=triangles_tri_super(i)
                ! print*,triangles_tri(nb_tri1)%p1%indice,triangles_tri(nb_tri1)%p2%indice,&
                ! triangles_tri(nb_tri1)%p3%indice,nb_tri1
            end if
        end do 

    call save_vtk(triangles_tri(1:nb_tri1),lespoints,"test.vtk",n+3,nb_tri1)


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Voronoi %%%%%%%%%
    allocate(lescentres(1:nb_tri1))
    do i=1,nb_tri1
        call construction_points_center(triangles_tri(i),lescentres(i))
    end do 


    nb_edg0=3*nb_tri1
    allocate(aretes_maille(1:nb_tri1,3),maille_aretes(1:nb_edg0,2),voronoi(nb_point,10))
    allocate(lesaretes(nb_edg0))
    allocate(maille_point(nb_point,10))
    maille_point=1
    call connectivite(lespoints,triangles_tri,lesaretes,aretes_maille,maille_aretes,maille_point,nb_point,nb_tri1,nb_aretes)
    ! do i=4,nb_point
    !     if (maille_point(i,1)<5) then
    !         print*,"point",i+3
    !         do j=2,maille_point(i,1)
    !             print*,triangles_tri(maille_point(i,j))%p1%indice,triangles_tri(maille_point(i,j))%p2%indice &
    !             ,triangles_tri(maille_point(i,j))%p3%indice
    !         end do 
    !     end if
    ! end do 
    ! do i=4,nb_point
    !     if (maille_point(i,1)<5) then
    !         voronoi(i,1)=maille_point(i,1)
    !         do j=2,maille_point(i,1)

    !             if(lesaretes(aretes_maille(nbr_tri_point,1))%p1%indice==lespoints(i)%indice)then 
    !                 aretes_deb=lesaretes(aretes_maille(nbr_tri_point,1))
    !             end if
    !             if(lesaretes(aretes_maille(nbr_tri_point,2))%p1%indice==lespoints(i)%indice)then 
    !                 aretes_deb=lesaretes(aretes_maille(nbr_tri_point,2))
    !             end if
    !             if(lesaretes(aretes_maille(nbr_tri_point,3))%p1%indice==lespoints(i)%indice)then 
    !                 aretes_deb=lesaretes(aretes_maille(nbr_tri_point,3))
    !             end if
    !         end do 
    !     end if

    ! end do 


    ! print*,"nombres des arrestes :", nb_aretes
    ! print*,"nombres des triangles :", nb_tri1

    ! do i=1,nb_aretes
    !     print*,"l'arret Numero :",i,"triangle à gauche",maille_aretes(i,1) ,"triangle à droit",maille_aretes(i,2) 
    ! end do

    ! print*,"******************",lesaretes(3)%p1%indice,lesaretes(3)%p2%indice

    ! do i=4,nb_point
    !         !do j=2,maille_point(i,1)
    !             indice_tri=maille_point(i,2)

    !             ! construire le new_tri
    !             call construire_new(triangles_tri(indice_tri),newtri,lespoints(i))
                
    !             !trouver les arretes de triangle indice_tri
    !             do j=1,3
    !                 edges(j)=lesaretes(aretes_maille(indice_tri,j))
    !             end do 
    !             !trouver un arrete de départ 
    !             call trouver_ar(triangles_tri(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
    !             do j=2,maille_point(i,1)



    !             end do 
                



    !         !end do 
    ! end do 


    do i=4,nb_point
            !do j=2,maille_point(i,1)
                indice_tri=maille_point(i,2)

                ! construire le new_tri
                call construire_new(triangles_tri(indice_tri),newtri,lespoints(i))
                
                !trouver les arretes de triangle indice_tri
                do j=1,3
                    edges(j)=lesaretes(aretes_maille(indice_tri,j))
                end do 
                !trouver un arrete de départ 
                call trouver_ar(triangles_tri(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
                do j=2,maille_point(i,1)



                end do 
                



            !end do 
    end do 
























    deallocate(lescentres)
    deallocate(aretes_maille,maille_aretes,maille_point,lesaretes)
    deallocate(lestriangles)
    deallocate(triangles_tri,triangles_tri_super)









contains
subroutine save_vtk (tri, points, name_file, n,nb_tri)
    type(triangle), dimension(:), intent(in)                :: tri
    type(point), dimension(:), intent(in)                   :: points            
    character(len=8), intent(in)                            :: name_file
    integer                                                 :: i
    integer, intent(in) :: n,nb_tri
    ! Open the file for writing
    open(unit=10, file=name_file, status='replace')
    ! Write informations
    write(10, '(A)') "# vtk DataFile Version 3.0"
    write(10, '(A)') "triangulation.vtk"
    write(10, '(A)') "ASCII"
    write(10, '(A)') "DATASET POLYDATA"
    ! Write vertex coordinates
    write(10, '(A,I0,A)') "POINTS ", n, " float"
    do i = 1, n
        write(10, '(2(F7.4,1X),F4.1)') points(i)%x , points(i)%y , 0.00000000
    end do
    ! Write cell connectivity
    write(10, '(A,I0,A,I0)') "POLYGONS ", nb_tri, " ", 4*nb_tri
    do i = 1, nb_tri
        write(10, '(4(I0,1X))') 3, tri(i)%p1%indice-1, tri(i)%p2%indice-1, tri(i)%p3%indice-1
    end do
    ! Close the file
    close(unit=10)
end subroutine

    
end program programme
