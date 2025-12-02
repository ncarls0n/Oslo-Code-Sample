! This programme is a quick check for 

program check_catalogue
    use :: iso_fortran_env, only : file_storage_size
    implicit none
    integer(kind=4) i, j, num_halos, num_halo_params, dummy
    real(kind=4)    R_THmax,z_thing
    real(kind=4),allocatable :: posxyz(:,:)  , velxyz(:,:) , &
                                posLxyz(:,:) , rth(:)      , F(:)
    character(len=128) catalogue,output
    character(len=:),allocatable :: tab
    tab='		'

    ! Read catalogue file name from command line
    read(*,*) catalogue

    ! Read catalogue file header and print it to the screen
    open(unit=1,file=catalogue,form='unformatted',access='stream')
    read(1) num_halos, R_THmax, z_thing
    close(1)
    write(*,*) 'Number of halos:  ' , num_halos
    write(*,*) 'Max filter scale: ' , R_THmax
    write(*,*) 'zthing:           ' , z_thing

    ! Allocate halo parameters:
    allocate(posxyz(3,num_halos))  ! Eulerian position
    allocate(velxyz(3,num_halos))  ! Velocity
    allocate(posLxyz(3,num_halos)) ! Lagrangian position
    allocate(rth(num_halos))       ! Filter scale
    allocate(F(num_halos))         ! Peak overdensity

    ! Read file size in bytes of halo catalogue, divide by 4 to get the
    ! number of 32-bit values in the catalogue, subtract the 3 values in
    ! the header, divide by the number of halos to get # 32-bit floats per
    ! halo
    inquire(file=catalogue, size=num_halo_params)
    num_halo_params=(num_halo_params/4-3)/num_halos
    write(*,*) 'The number of halo parameters is ', num_halo_params

    ! Read catalogue main body
    open(unit=1,file=catalogue,form='unformatted',access='stream')
    read(1) dummy, dummy, dummy
    read(1) ( ( posxyz(j,i),  j=1,3                ),         &
              ( velxyz(j,i),  j=1,3                ), rth(i), &
              ( posLxyz(j,i), j=1,3                ), F(i),   &
              ( dummy,        j=12,num_halo_params ), i=1,num_halos )
    close(1)

    write(*,*) 'x'//tab//'y'//tab//'z'//tab//'v_x'//tab//'v_y'//tab//'v_z'&
        //tab//'R_TH'//tab//'x_L'//tab//'y_L'//tab//'z_L'//tab//'F_pk'
    
    write(*,*) ( ( posxyz(j,i),  j=1,3 ),         &
                 ( velxyz(j,i),  j=1,3 ), rth(i), &
                 ( posLxyz(j,i), j=1,3 ), F(i),   i=1,2 ) 
    write(*,*) '.'
    write(*,*) '.'
    write(*,*) '.'
    write(*,*) ( ( posxyz(j,i),  j=1,3 ),         &
                 ( velxyz(j,i),  j=1,3 ), rth(i), &
                 ( posLxyz(j,i), j=1,3 ), F(i),   i=num_halos-1,num_halos )

end program check_catalogue
