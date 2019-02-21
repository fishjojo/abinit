#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_paw_mkvij

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_xmpi

 use m_pawrhoij,   only : pawrhoij_type,pawrhoij_init_unpacked
 use m_pawcprj,    only : pawcprj_type,pawcprj_alloc,pawcprj_get,pawcprj_free
 use m_paral_atom, only : get_my_atmtab,free_my_atmtab


contains


!!****f* m_paw_mkvij/pawmkvij
!!
!! NAME
!! pawmkvij
!!
!! FUNCTION
!! Compute vemb related term; similar to pawmkrhoij
!!
!! TODO
!! Parallelization
!!
!! SOURCE

subroutine pawmkvij(atindx,atindx1,vemb,norb,cprj,dimcprj,istwfk,mband,mband_cprj,mcprj,mkmem,mpi_enreg,&
&                      natom,nband,nkpt,nspinor,nsppol,paral_kgb,&
&                      pawrhoij,unpaw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkvij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: norb
 integer,intent(in) :: mband,mband_cprj,mcprj,mkmem,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: paral_kgb,unpaw
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom),istwfk(nkpt)
 integer,intent(in) :: nband(nkpt*nsppol)
 type(pawcprj_type),target,intent(in) :: cprj(natom,mcprj)
 type(pawrhoij_type),intent(inout),target:: pawrhoij(:)   !simply use pawrhoij as pawvij for now
 complex(dpc),intent(in) :: vemb(norb,norb) 

!Local variables ---------------------------------------
!scalars
 integer :: ibg,ikpt,iorder_cprj,isppol,nband_k,nrhoij
 integer :: me,spaceComm,my_nspinor,cplex
 logical :: paral_atom
 character(len=500) :: msg

!arrays
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 type(pawcprj_type),pointer :: cprj_ptr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_all(:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

!Init MPI data
 spaceComm=mpi_enreg%comm_kpt
 me=mpi_enreg%me_kpt

!Check size of cprj
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 if (mcprj/=my_nspinor*mband_cprj*mkmem*nsppol) then
   msg=' wrong size for cprj !'
   MSG_BUG(msg)
 end if

!Size of pawrhoij datastructure
 nrhoij=size(pawrhoij)

!Check if pawrhoij is distributed over atomic sites
 paral_atom=(nrhoij/=natom.and.mpi_enreg%nproc_atom>1)
 if (paral_atom.and.nrhoij/=mpi_enreg%my_natom) then
   msg=' Size of pawrhoij should be natom or my_natom !'
   MSG_BUG(msg)
 end if

!Allocate temporary cwaveprj storage
 if(nspinor.gt.1) MSG_BUG(' nspinor > 1 NYI!')
 ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,nspinor*norb))
 call pawcprj_alloc(cwaveprj,0,dimcprj)

!Initialize temporary file (if used)
 iorder_cprj=0

!Build and initialize unpacked rhoij (to be computed here)
 call pawrhoij_init_unpacked(pawrhoij)

!If pawrhoij is MPI-distributed over atomic sites, gather it
 if (paral_atom) then
   ABI_DATATYPE_ALLOCATE(pawrhoij_all,(natom))
 else
   pawrhoij_all => pawrhoij
 end if

 if(nkpt > 1) MSG_BUG(' nkpt > 1 NYI!')

!LOOP OVER SPINS
 ibg=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   do ikpt=1,nkpt

     cplex = 2;if (istwfk(ikpt)>1) cplex=1
     nband_k=nband(ikpt+(isppol-1)*nkpt)

     cprj_ptr => cprj

     call pawcprj_get(atindx1,cwaveprj,cprj_ptr,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&      mband_cprj,mkmem,natom,norb,nband_k,nspinor,nsppol,unpaw,&
&      mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

     call pawaccvij(atindx,cplex,vemb,norb,cwaveprj,cwaveprj,isppol,nrhoij,natom,nspinor,pawrhoij_all)

   enddo !ikpt
 enddo !isppol

 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

 DBG_EXIT("COLL")

end subroutine pawmkvij



!!****f* m_paw_mkvij/pawaccvij
!!
!! NAME
!! pawaccvij
!!
!! FUNCTION
!! Compute vemb related term; similar to pawaccrhoij
!!
!! TODO
!! Parallelization
!!
!! SOURCE
subroutine pawaccvij(atindx,cplex,vemb,norb,cwaveprj,cwaveprj1,isppol,my_natom,natom,&
&                       nspinor,pawrhoij,&
&                       comm_atom,mpi_atmtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawaccvij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: norb
 integer,intent(in) :: cplex,isppol,my_natom,natom,nspinor
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,intent(in) :: atindx(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawcprj_type),intent(in) :: cwaveprj(natom,nspinor),cwaveprj1(natom,nspinor)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom) !simply use pawrhoij as pawvij for now
 complex(dpc),intent(in) :: vemb(norb,norb)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iatm,iatom,iatom1,ilmn,iplex,j0lmn,jlmn,klmn,iorb,jorb
 integer :: my_comm_atom
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: ro11_re
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: cpi0(2,nspinor*norb),cpj0(2,nspinor*norb),cpjv(2,nspinor*norb)
 real(dp) :: vij_re(norb,norb),vij_im(norb,norb)

! ***********************************************************************

 DBG_ENTER("COLL")

 if (my_natom==0) return

 vij_re = real(vemb,kind=dp)
 vij_im = aimag(vemb)


!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 if(paral_atom) MSG_BUG(' paral_atom is not implemented yet!')
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

 if(cplex/=2) MSG_BUG(' cwaveprj should be complex!')

 if (nspinor==1) then
     cplex_rhoij=pawrhoij(1)%cplex
     if (cplex_rhoij==1) then
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2

           do iorb=1,norb
             cpj0(1:cplex,iorb)=cwaveprj(iatm,iorb)%cp(1:cplex,jlmn)
           enddo
           
           do jorb=1,norb 
             cpjv(1,jorb) = dot_product(cpj0(1,:), vij_re(:,jorb) ) - dot_product(cpj0(2,:), vij_im(:,jorb) )
             cpjv(2,jorb) = dot_product(cpj0(1,:), vij_im(:,jorb) ) + dot_product(cpj0(2,:), vij_re(:,jorb) )
           enddo

           do ilmn=1,jlmn
             klmn=j0lmn+ilmn

             do jorb=1,norb
               cpi0(1:cplex,jorb)=cwaveprj1(iatm,jorb)%cp(1:cplex,ilmn)
             enddo

             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dot_product( cpjv(iplex,:),cpi0(iplex,:) )
             end do

             pawrhoij(iatom)%rhoij_(klmn,isppol)=pawrhoij(iatom)%rhoij_(klmn,isppol)+ro11_re
           end do
         end do
       end do
     else
       MSG_BUG(' pawvij has to be real!')
     end if
 else !nspinor==2
     MSG_BUG(' nspinor > 1 is not implemented yet!')
 endif

 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")


end subroutine pawaccvij


end module m_paw_mkvij
