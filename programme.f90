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
    type(point),allocatable,dimension(:)::lespoints,lescentres,point_delaunay,voronoi_entre,voronoi_sort,point_ajouter
    type(edge),allocatable,dimension(:)::lesaretes,aretes_delaunay
    type(triangle),allocatable,dimension(:)::lestriangles,triangles_tri_super,triangles_tri,tri_delaunay
    integer,allocatable, dimension(:,:):: aretes_maille,maille_aretes,voronoi
    integer,allocatable,dimension(:,:)::maille_point,centre_bord
    integer,allocatable,dimension(:)::point_bord1,edge_bord1,point_bord,edge_bord,poin_tri,poin_tri1
    type(edge)::aretes_deb1,arretes_ari,aretes_deb
    type(point)::p1,p2,p3
    type(edge)::edge1,edge2,edge3
    type(triangle)::newtri,newtri1
    type(edge),dimension(3)::edges
    integer::indice_ag,jj,nb_cel,countour,indice_tri1, compteur_point, compteur_edge,compteur_pointtri,nb_point_aj
    logical::bord


    
    !nombre des points pour chaque intervalle 
    np=10
    n=np*np
    nb_point=n+3
    ! Bornes de l'intervalle
    xmin = -5.0
    xmax = 5.0
    ymin = -5.0
    ymax = 5.0

    allocate(lespoints(400000))
    !initialiser la nuage de point
    !call nuage(xmax, xmin, ymax, ymin, np, lespoints)
    call nuage_trigo_v2(nb_point,lespoints,ymax,37)
   ! print*,"c'est bon"
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
                ! print*,triangles_tri(nb_tri1)%p1%indice,triangles_tri(nb_tri1)%p2%indice,&triangles_tri(nb_tri1)%p3%indice,nb_tri1
            end if
        end do 

    call save_vtk(triangles_tri(1:nb_tri1),lespoints,"test.vtk",nb_point,nb_tri1)
!     allocate(point_ajouter(nb_point))
!     nb_point_aj=0
!         do i=1,nb_tri1
!             if(hasObtuseAngle(triangles_tri(i)%p1,triangles_tri(i)%p2,triangles_tri(i)%p3) .eqv. .true.) then
!                 nb_point_aj=nb_point_aj+1
!                 point_ajouter(nb_point_aj)=getMidpointOfObtuseSide(triangles_tri(i)%p1,triangles_tri(i)%p2,triangles_tri(i)%p3)
!                 !print*,"errrrrrrrrr",triangles_tri(i)%p1%indice,triangles_tri(i)%p2%indice,triangles_tri(i)%p3%indice
!             end if 
!         end do 
!     print*,"point ajouter ",nb_point_aj
!    ! call removeDuplicates(poin_tri(1:compteur_pointtri),poin_tri1,compteur_pointtri)
!     allocate(point_delaunay(nb_point+nb_point_aj))

!     do i=1,nb_point
!         !print*,"errrrrrr",poin_tri1(i),i,compteur_pointtri
!         point_delaunay(i)=lespoints(i)
!         point_delaunay(i)%indice=i
!     end do 
!     do i=1,nb_point_aj
!             point_delaunay(i+nb_point)=point_ajouter(i)
!             point_delaunay(i+nb_point)%indice=i+nb_point
!     end do
!     deallocate(lespoints)
!     nb_point=nb_point+nb_point_aj
!     allocate(lespoints(nb_point))
!     lespoints=point_delaunay

!     deallocate(point_delaunay,triangles_tri,lestriangles,triangles_tri_super,point_ajouter)

! !%%%%%%%%%%%%%%%%%%%%%%%hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
!     super = SuperTriangle(lespoints)
!     call construction_points(lespoints(1),super%p1%x,super%p1%y,1)
!     call construction_points(lespoints(2),super%p2%x,super%p2%y,2)
!     call construction_points(lespoints(3),super%p3%x,super%p3%y,3)
    
!     allocate(lestriangles(3*nb_point))

!     do i=1,3*nb_point
!         lestriangles(i)=super
!     end do

!     call algo_tri(lespoints,lestriangles,nb_point)

!     !trier les triangles 
!     allocate(triangles_tri_super(3*nb_point))
!     nb_tri=1
!     triangles_tri_super(1)=super
!     do i=1,3*nb_point
!        if(lestriangles(i)%p1%indice/=1 .or. lestriangles(i)%p2%indice/=2 &
!        .or. lestriangles(i)%p3%indice/=3)then
!         triangles_tri_super(i+1)=lestriangles(i)
!             nb_tri=nb_tri+1
!        end if 
!     end do 
!     !call save_vtk(triangles_tri_super(1:nb_tri),lespoints,"testt.vtk",n+3,nb_tri)

!     allocate(triangles_tri(3*nb_tri))
!         nb_tri1=0
!         do i=1,nb_tri
!             if(triangles_tri_super(i)%p1%indice/=1 .and. triangles_tri_super(i)%p1%indice/=2 &
!             .and. triangles_tri_super(i)%p1%indice/=3 .and. triangles_tri_super(i)%p2%indice/=1 &
!             .and. triangles_tri_super(i)%p2%indice/=2 .and. triangles_tri_super(i)%p2%indice/=3 &
!             .and. triangles_tri_super(i)%p3%indice/=1 .and. triangles_tri_super(i)%p3%indice/=2 &
!             .and. triangles_tri_super(i)%p3%indice/=3)then
!                 nb_tri1=nb_tri1+1
!                 triangles_tri(nb_tri1)=triangles_tri_super(i)
!                 ! print*,triangles_tri(nb_tri1)%p1%indice,triangles_tri(nb_tri1)%p2%indice,&
!                 ! triangles_tri(nb_tri1)%p3%indice,nb_tri1
!             end if
!         end do 

!     call save_vtk(triangles_tri(1:nb_tri1),lespoints,"tett.vtk",nb_point,nb_tri1)



!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% netoyage
allocate(tri_delaunay(nb_tri1))
    do i=1,nb_tri1
        tri_delaunay(i)=triangles_tri(i)
        tri_delaunay(i)%indice=i
    end do 
deallocate(lestriangles)
deallocate(triangles_tri,triangles_tri_super)

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Voronoi %%%%%%%%%


    nb_edg0=3*nb_tri1
    allocate(aretes_maille(1:nb_tri1,3),maille_aretes(1:nb_edg0,2),voronoi(nb_point,10))
    allocate(lesaretes(nb_edg0))
    allocate(maille_point(nb_point,20))
    maille_point=0
    call connectivite(lespoints,tri_delaunay,lesaretes,aretes_maille,maille_aretes,maille_point,nb_point,nb_tri1,nb_aretes)
    
    open(unit=10, file="mail.txt", status='replace')
    do i=1,nb_point
        do j=1,10
            write(10, '(4(I0,1X))', ADVANCE='NO') maille_point(i,j)
        end do 
        write(10,*)
    end do 

    close(10)

    allocate(aretes_delaunay(nb_aretes),edge_bord1(nb_aretes),point_bord1(nb_aretes))
        !trouver les arrtes et les points au bord
    compteur_edge = 0
    compteur_point=0
    do i=1,nb_aretes
        if (maille_aretes(i,1) == 0 .or. maille_aretes(i,2) == 0) then 
            compteur_edge = compteur_edge + 1
            edge_bord1(compteur_edge)=i
            compteur_point=compteur_point+1
            point_bord1(compteur_point)=lesaretes(i)%p1%indice
            compteur_point=compteur_point+1
            point_bord1(compteur_point)=lesaretes(i)%p2%indice
        end if 
    end do 

    call removeDuplicates(point_bord1(1:compteur_point),point_bord,compteur_point)
    allocate(edge_bord(1:compteur_edge))
    edge_bord=edge_bord1(1:compteur_edge)

    deallocate(point_bord1,edge_bord1)

    print *,"points au bord",size(point_bord)
    print *,"edge au bord",compteur_edge
    
    allocate(centre_bord(size(edge_bord),4))
    do i=1,size(edge_bord)
        centre_bord(i,1)=edge_bord(i)
        centre_bord(i,2)=lesaretes(edge_bord(i))%p1%indice
        centre_bord(i,3)=lesaretes(edge_bord(i))%p2%indice
    end do

    allocate(lescentres(1:nb_tri1+size(point_bord)))
    do i=1,nb_tri1
        call construction_points_center(tri_delaunay(i),lescentres(i))
    end do 
    do i=nb_tri1+1,nb_tri1+size(edge_bord)
        lescentres(i)=findMidpoint(lespoints(centre_bord(i-nb_tri1,2)),lespoints(centre_bord(i-nb_tri1,3)))
        lescentres(i)%indice=i
        centre_bord(i-nb_tri1,4)=i
    end do 


    do i=1,nb_aretes
        aretes_delaunay(i)=lesaretes(i)
    end do 
    deallocate(lesaretes)
  !  print *, "je suis ici"
    
    voronoi=0
    do i=4,nb_point
        bord=.false.
        do j=1,size(point_bord)
            if(i==point_bord(j)) then
                bord=.true.
            end if 
        end do 
        if(bord .eqv. .false. ) then 
            voronoi(i,1)=maille_point(i,1)
            !print*,"point ",i
            allocate(voronoi_entre(maille_point(i,1)),voronoi_sort(maille_point(i,1)))
            do j=1,maille_point(i,1)
                voronoi_entre(j)=lescentres(maille_point(i,j+1))
                !write(*,'(i3)', advance='no')voronoi_entre(j)%indice
            end do
            !write(*,*) "" 
            call sortPointsCCW(lespoints(i),voronoi_entre,voronoi_sort)
            !print*,"point ",i
            write(*,'(i3)', advance='no')maille_point(i,1)
            do j=1,size(voronoi_entre)
              !  write(*,'(i3)', advance='no')voronoi_sort(j)%indice-1
                voronoi(i,j+1)=voronoi_sort(j)%indice-1
            end do 
           ! write(*,*) "" 
            deallocate(voronoi_entre,voronoi_sort)
        end if 
        if(bord .eqv. .true.) then
            voronoi(i,1)=maille_point(i,1)+2
            allocate(voronoi_entre(maille_point(i,1)+2),voronoi_sort(maille_point(i,1)+2))
            do j=1,maille_point(i,1)
                voronoi_entre(j)=lescentres(maille_point(i,j+1))
                !write(*,'(i3)', advance='no')voronoi_entre(j)%indice
            end do
            compteur_point=0
            do j=1,size(edge_bord)
                if(centre_bord(j,2)==i .or. centre_bord(j,3)==i  ) then 
                    compteur_point=compteur_point+1
                    voronoi_entre(maille_point(i,1)+compteur_point)=lescentres(centre_bord(j,4))
                end if
                if(compteur_point>2) then 
                    print*,"j'ai un probleme au bord ligne 275"
                end if 
            end do 
          !  write(*,*) "" 
            call sortPointsCCW(lespoints(i),voronoi_entre,voronoi_sort)
            !print*,"point ",i
           ! write(*,'(i3)', advance='no')maille_point(i,1)
            do j=1,size(voronoi_entre)
               ! write(*,'(i3)', advance='no')voronoi_sort(j)%indice-1
                voronoi(i,j+1)=voronoi_sort(j)%indice-1
            end do 
          !  write(*,*) "" 
            deallocate(voronoi_entre,voronoi_sort)

        end if 


    end do 
 
    


    ! voronoi=0
    ! do i=4,nb_point
    !     countour=0
    !     !verifiez si le point est au bord 
    !     do j=2,maille_point(i,1)
    !         ! print *, "je suis ici"
    !         indice_tri=maille_point(i,j)
    !         do jj=1,3
    !             if(maille_aretes(aretes_maille(indice_tri,jj),1)==0) then
    !                 ! if (lesaretes(aretes_maille(indice_tri,jj))%p1%indice==lespoints(i)%indice .or. &
    !                 !   lesaretes(aretes_maille(indice_tri,jj))%p2%indice==lespoints(i)%indice) then 
    !                 !         countour=countour+1
    !                 ! end if
    !                 print*,"888888888888",lesaretes(aretes_maille(indice_tri,jj))%p1%indice
    !             end if
    !             if(maille_aretes(aretes_maille(indice_tri,jj),2)==0) then
    !                 if (lesaretes(aretes_maille(indice_tri,jj))%p1%indice==lespoints(i)%indice .or. &
    !                 lesaretes(aretes_maille(indice_tri,jj))%p2%indice==lespoints(i)%indice) then 
    !                     countour=countour+1
    !                 end if
    !             end if
    !         end do
    !     end do 

    !     print*,"***********************",countour
    !     if( countour/=2)then
    !             print*, "le point ",i
                
    !             do j=2,maille_point(i,1)+1
    !                 indice_tri=maille_point(i,j)
    !                 print*,"maille",indice_tri
    !                 call construire_new(tri_delaunay(indice_tri),newtri,lespoints(i))
    !                 do jj=2,maille_point(i,1)+1
    !                     if(jj/=j) then
    !                         indice_tri1=maille_point(i,jj)
    !                         call construire_new(tri_delaunay(indice_tri1),newtri1,lespoints(i))
    !                         if(newtri1%p2%indice==newtri%p3%indice ) then
    !                             print *,"maille voisine",indice_tri, indice_tri1
    !                         end if
    !                     end if
    !                 end do 



    !             end do 

        

    !     end if 



        






    ! end do 



!     voronoi=0
!     do i=4,nb_point
!                 indice_tri=maille_point(i,2)
!                 print*,"premeir indice ",indice_tri,tri_delaunay(maille_point(i,2))%indice

!                 voronoi(i,1)=maille_point(i,1)
!                 !trouver les arretes de triangle indice_tri
!                 print*,indice_tri
!                 do j=1,3
!                     edges(j)=aretes_delaunay(aretes_maille(indice_tri,j))
!                     print*,edges(j)%p2%indice,edges(j)%p1%indice,lespoints(i)%indice
!                 end do 
                
! !!!! c'est bon 
!                 !trouver un arrete de départ 
!                 call trouver_ar(tri_delaunay(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
!                 print*,"5555555555555555555555555555555555",aretes_deb%indice,arretes_ari%indice
!                 do j=2,maille_point(i,1)
!                     indice_tri=maille_point(i,j)
!                     indice_ag=maille_aretes(aretes_deb%indice,1)
!                     voronoi(i,j)=indice_ag

!                     if(j/=maille_aretes(i,1)) then 
!                         indice_tri=maille_point(i,j+1)
!                         !trouver les arretes de triangle indice_tri
!                         do jj=1,3
!                             edges(jj)=aretes_delaunay(aretes_maille(indice_tri,jj))
!                         end do 
!                         !trouver un arrete de départ 
!                         call trouver_ar(tri_delaunay(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
!                     end if 
!                 end do 
!    end do 
   nb_cel=0
   do i=4,nb_point
   nb_cel=nb_cel+voronoi(i,1)
   end do 
   nb_cel=nb_cel+nb_point-3

   open(unit=10, file="vorn.vtk", status='replace')
   ! Write informations
   write(10, '(A)') "# vtk DataFile Version 3.0"
   write(10, '(A)') "triangulation.vtk"
   write(10, '(A)') "ASCII"
   write(10, '(A)') "DATASET POLYDATA"
   ! Write vertex coordinates
   write(10, '(A,I0,A)') "POINTS ", nb_tri1+size(edge_bord), " float"
   do i = 1, nb_tri1+size(edge_bord)
       write(10, '(2(F10.4,1X),F4.1)') lescentres(i)%x , lescentres(i)%y , 0.00000000
   end do
   ! Write cell connectivity
   write(10, '(A,I0,A,I0)') "POLYGONS ", nb_point-3, " ", nb_cel
   do i = 4,nb_point
        do j=1,voronoi(i,1)+1
            write(10, '(4(I0,1X))', ADVANCE='NO') voronoi(i,j)
        end do 
        write(10,*)
   end do
   ! Close the file
   close(unit=10)

   !##############"""""""""groupe CDO
   open(unit=10, file="meshD.txt", status='replace')    
   write(10,*) nb_point-3
   print*,nb_point-3
   do i=4,nb_point 
       write(10,*) lespoints(i)%x,lespoints(i)%y
       print*,lespoints(i)%x,lespoints(i)%y
   end do

   write(10,*) nb_tri1
   do i=1,nb_tri1
       write(10,'(1(I0,1X))') 3
       write(10,'(3(I0,1X))') tri_delaunay(i)%p1%indice-3,tri_delaunay(i)%p2%indice-3&
       ,tri_delaunay(i)%p3%indice-3
   end do 
   close(10)

   open(unit=10, file="meshv.txt", status='replace')    
   write(10,*) nb_tri1+size(edge_bord)
   do i=1,nb_tri1+size(edge_bord)
       write(10,*) lescentres(i)%x,lescentres(i)%y
   end do

   write(10,*) nb_point-3
    do i = 4,nb_point
        write(10, '(4(I0,1X))') voronoi(i,1)
        do j=2,voronoi(i,1)+1
            write(10, '(4(I0,1X))', ADVANCE='NO') voronoi(i,j)+1
        end do 
        write(10,*)
   end do
   ! Close the file
   close(unit=10)

   !!ecrire maillage vornio





    ! do i=1,nb_aretes
    !     print*,"l'arret Numero :",i,"triangle à gauche",maille_aretes(i,1) ,"triangle à droit",maille_aretes(i,2) 
    ! end do

    ! print*,"******************",aretes_delaunay(3)%p1%indice,aretes_delaunay(3)%p2%indice

    ! do i=4,nb_point
    !         !do j=2,maille_point(i,1)
    !             indice_tri=maille_point(i,2)

    !             ! construire le new_tri
    !             call construire_new(tri_delaunay(indice_tri),newtri,lespoints(i))
                
    !             !trouver les arretes de triangle indice_tri
    !             do j=1,3
    !                 edges(j)=lesaretes(aretes_maille(indice_tri,j))
    !             end do 
    !             !trouver un arrete de départ 
    !             call trouver_ar(tri_delaunay(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
    !             do j=2,maille_point(i,1)



    !             end do 
                



    !         !end do 
    ! end do 



!     do i=4,6
!             !do j=2,maille_point(i,1)
!                 indice_tri=tri_delaunay(maille_point(i,2))%indice
!                 voronoi(i,1)=maille_point(i,1)-1
!                 !trouver les arretes de triangle indice_tri
!                 print*,indice_tri
!                 do j=1,3
!                     edges(j)=aretes_delaunay(aretes_maille(indice_tri,j))
!                     print*,edges(j)%p2%indice,edges(j)%p1%indice,lespoints(i)%indice
!                 end do 
                

!                 !trouver un arrete de départ 
!                 ! call trouver_ar(tri_delaunay(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
!                 ! aretes_deb=aretes_deb1
!                 ! print*,"5555555555555555555555555555555555",aretes_deb%indice
!                 ! do j=2,maille_point(i,1)
!                 !     indice_tri=maille_point(i,j)
!                 !     !indice_ag=maille_aretes(aretes_deb%indice,1)
!                 !     voronoi(i,j)=indice_ag

!                 !     if(j/=maille_aretes(i,1)) then 
!                 !         indice_tri=maille_point(i,j+1)
!                 !         !trouver les arretes de triangle indice_tri
!                 !         do jj=1,3
!                 !             edges(jj)=aretes_delaunay(aretes_maille(indice_tri,jj))
!                 !         end do 
!                 !         !trouver un arrete de départ 
!                 !         call trouver_ar(tri_delaunay(indice_tri),edges,lespoints(i),aretes_deb,arretes_ari)
!                 !     end if 
!                 ! end do 
                



!     !         !end do 
!    end do 
!     do i=4,nb_point
!         print *, "le point d'indice",i
!         print *,"nb de maille ", maille_point(i,1)-1
!         print *,"1ere maille", maille_point(i,2)
!         do j=2,maille_point(i,1)
!             print*,"****************maille d'indice ",maille_point(i,j)
!             print*,tri_delaunay(maille_point(i,j))%p1%indice,tri_delaunay(maille_point(i,j))%p2%indice&
!             ,tri_delaunay(maille_point(i,j))%p3%indice
!             print*,"etape2"
!             print*,"edge1",aretes_delaunay(aretes_maille(maille_point(i,j),1))%p1%indice,&
!             aretes_delaunay(aretes_maille(maille_point(i,j),1))%p2%indice
!             print*,"edge2",aretes_delaunay(aretes_maille(maille_point(i,j),2))%p1%indice,&
!             aretes_delaunay(aretes_maille(maille_point(i,j),2))%p2%indice
!             print*,"edge3",aretes_delaunay(aretes_maille(maille_point(i,j),3))%p1%indice,&
!             aretes_delaunay(aretes_maille(maille_point(i,j),3))%p2%indice
        
!         print*, "'''''''''''''''''''''''''''"
!         indice_tri=maille_point(i,j)
!         voronoi(i,1)=maille_point(i,1)-1
!         !trouver les arretes de triangle indice_tri
!         print*,indice_tri
!             do jj=1,3
!                 edges(jj)=aretes_delaunay(aretes_maille(indice_tri,jj))
!                 print*,edges(jj)%p2%indice,edges(jj)%p1%indice,lespoints(i)%indice
!             end do 

!          end do 



!     end do 


!     print*,"#################################################"

!     do i=1,nb_tri1
!         print*,i,tri_delaunay(i)%p1%indice
!         print*,i,tri_delaunay(i)%indice




!     end do 




















    deallocate(lescentres)
    deallocate(aretes_maille,maille_aretes,maille_point,aretes_delaunay,tri_delaunay)










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
