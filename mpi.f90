
subroutine  update_mpi_pos ! send the position from proc 0 to all others
use mpi
use system
implicit none
CALL MPI_BCAST(xpos, di*Npart, MPI_REAL, 0, MPI_COMM_WORLD,err)
end

subroutine update_mpi_force ! send the forces from all procs to proce 0
use mpi
use BD
use system
implicit none
forces = 0.0
CALL MPI_REDUCE(forces_tosend, forces, di*Npart, MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, err)
end

subroutine update_mpi_pres ! send the forces from all procs to proce 0
use mpi
use BD
use system
implicit none
pres = 0.0
CALL MPI_ALLREDUCE(pres_tosend, pres, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, err)
end

subroutine update_mpi_energy ! send the forces from all procs to proce 0
use mpi
use BD
use system
implicit none
senergy = 0.0
CALL MPI_REDUCE(senergy_tosend, senergy, 1, MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, err)
end

subroutine proc
! decompose cells into processors for parallelization
use system
use nlist
use mpi
implicit none

integer j, i

if(size.gt.ncell) then
print*, 'More processors than cells...'
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif 

allocate(Nlistproc(size))
allocate(listproc(size, int(ncell/size)+1))

i = 1
Nlistproc(:) = 0
listproc(:,:) = 0

do j = 1, ncell
Nlistproc(i) = Nlistproc(i)+1 
listproc(i, Nlistproc(i)) = j
i = i + 1 
if(i.gt.size)i=1
enddo

if(rank.eq.0)print*, 'Proc', size,'cells/Proc', int(ncell/size)+1

!do i = 1, size
!if(rank.eq.0)print*, 'Processor ',  i, ' does cells'
!do j = 1, Nlistproc(i) 
!if(rank.eq.0)print*, listproc(i, j)
!enddo
!enddo

end


subroutine update_mpi_bo ! send bo from all processors
use mpi
use BD
use system
use bo
use grsystem
implicit none
CALL MPI_REDUCE(bo4r_tosend, bo4r, types*types*Ngrmax, MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(bo2r_tosend, bo2r, types*types*Ngrmax, MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(bo6r_tosend, bo6r, types*types*Ngrmax, MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(bo4rt_tosend, bo4rt, Ngrmax, MPI_REAL, MPI_SUM, MPI_COMM_WORLD,0, err)
CALL MPI_REDUCE(bo2rt_tosend, bo2rt, Ngrmax, MPI_REAL, MPI_SUM, MPI_COMM_WORLD,0, err)
CALL MPI_REDUCE(bo6rt_tosend, bo6rt, Ngrmax, MPI_REAL, MPI_SUM, MPI_COMM_WORLD,0, err)
end

subroutine update_mpi_gr ! send bo from all processors
use mpi
use BD
use system
use grsystem
implicit none
CALL MPI_REDUCE(GR_tosend, GR, types*types*Ngrmax, MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, err)
end

subroutine update_mpi_bo_vecinos ! send bo from all processors
use mpi
use BD
use system
use bo
use grsystem
implicit none
CALL MPI_REDUCE(NVECINOS2_tosend, NVECINOS2, Npart, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(VECINOS2_tosend, VECINOS2, Npart*VV, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(NVECINOS4_tosend, NVECINOS4, Npart, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(VECINOS4_tosend, VECINOS4, Npart*VV, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(NVECINOS6_tosend, NVECINOS6, Npart, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(VECINOS6_tosend, VECINOS6, Npart*VV, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)
CALL MPI_REDUCE(NVECINOS1_tosend, NVECINOS1, Npart, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD, err)
CALL MPI_REDUCE(VECINOS1_tosend, VECINOS1, Npart*VV, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD, err)
CALL MPI_REDUCE(NVECINOS3_tosend, NVECINOS3, Npart, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD, err)
CALL MPI_REDUCE(VECINOS3_tosend, VECINOS3, Npart*VV, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD, err)


end

subroutine update_mpi_bo_aux ! send bo from all processors
use mpi
use BD
use system
use bo
use grsystem
implicit none
CALL MPI_ALLREDUCE(bo2aux_tosend, bo2aux, Npart, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, err)
CALL MPI_ALLREDUCE(bo4aux_tosend, bo4aux, Npart, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, err)
CALL MPI_ALLREDUCE(bo6aux_tosend, bo6aux, Npart, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, err)
CALL MPI_ALLREDUCE(bo1aux_tosend, bo1aux, Npart, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, err)
CALL MPI_ALLREDUCE(bo3aux_tosend, bo3aux, Npart, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, err)
end


