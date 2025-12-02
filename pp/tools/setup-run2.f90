! setup-run.f90
! Created by: Nathan J. Carlson
! Last edit:  29 November 2022
! 
! USAGE:
! 

program setuprun
    use setuprunfuncs
    implicit none
    integer n, iLPT, cat_len, n_proc, totalnodes, n_buff_min, n_buff, n_tile, n_mesh_max, n_fields, n_eff
    real    s, s_buff, noderam, s_latt
    character(len=7) cluster

    read(*,*) n
    read(*,*) s
    read(*,*) cluster
    read(*,*) s_buff
    read(*,*) iLPT
    read(*,*) cat_len

    if(cluster/='Niagara') then
        write(*,*) 'Peak Patch currently only runs on Niagara'
        stop
    elseif(cluster=='Niagara') then
        n_proc      = 40   ! number of processors/node on cluster
        noderam    = 188. ! RAM/node in GiB
        totalnodes = 1000 ! available nodes on cluster
    endif

    if(iLPT==1) then
        n_fields=4
    elseif(iLPT==2) then
        n_fields=7
    else
        write(*,*) 'Only setup to do 1LPT or 2LPT'
        stop
    endif

    n_buff_min = int( s_buff/(s+2*s_buff)*float(n) )+1

    n_mesh_min = 128
    n_mesh_max = int(     &! 
        (0.5*noderam      &! 1/2 RAM in GiB per node
         /n_proc          &! / number of CPUs per node
         *2**30           &! * 2^30 B/GiB
         /4               &! / 4 B per 32-b float
         /(n_fields+1.25) &! / ( 7 fields + 1 working field + 1/4 extra )
        )**(1./3.)        &! cube root because n_mesh^3 is # of voxels
        +.5)              &! round to the nearest integer

    do n_tile=1,n/n_buff_min-2
        n_buff = n_buff_min
        n_eff  = get_n_eff(n,n_buff)        
        s_latt = get_s_latt(n_eff,s)

        

    enddo

! Need to find:
! n_tile
! n_mesh
! n_buff
! nnodes
! tpnode min(ntile^3,n_proc)

! Conditions that must be met:
! mod(n,n_tasks)==0 where n_tasks = n_nodes*n_tasks_per_node

end program setuprun

module setuprunfuncs
    implicit none
    contains

        integer function get_n_eff(n_in, n_buff_in)
            integer, intent(in) :: n_in, n_buff_in
            get_n_eff=n_in-2*n_buff_in
        end function get_n_eff

        real function get_s_latt(n_eff_in, s_in)
            integer, intent(in) :: n_eff_in
            real, intent(in)    :: s_in
            get_s_latt=s_in/float(n_eff_in)
        end function get_s_latt

        integer function fftw_mem(n_in,n_nodes_in,n_fields_in,n_nodes_in)
            integer, intent(in) :: n_in,n_nodes_in,n_fields_in,n_nodes_in
            fftw_mem=n_fields    &! 7 fields per tile
                     *4          &! * number of B per 32-b float
                     /2**30      &! / number of B per GiB
                     *n**3       &! * n^3
                     *n_nodes_in &! * number of nodes

!        function get_n_buff
!            
!        end function get_n_buff

end module setuprunfuncs
