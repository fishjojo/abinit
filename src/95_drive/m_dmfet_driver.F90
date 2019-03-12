!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dmfet_driver
!! NAME
!!  m_dmfet_driver
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2018- Xing Zhang
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_dmfet_driver

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_dmfet
 use m_hdr
 use m_ebands
 use m_errors
 use m_xmpi
 use m_bandfft_kpt

 use m_symtk,            only : matr3inv
 use m_common,           only : setup1
 use m_kg,               only : kpgio
 use m_crystal,          only : crystal_t, crystal_free, crystal_print
 use m_crystal_io,       only : crystal_from_hdr
 use defs_wvltypes,      only : wvl_data
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawfgr,           only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_initylmg,         only : initylmg
 use m_pspini,           only : pspini
 use m_wfk,              only : wfk_t,wfk_read_eigenvalues,wfk_open_read,wfk_read_band_block
 use m_io_tools,         only : iomode_from_fname,get_unit

 implicit none

 private
!!***

 public :: dmfet
!!***

contains
!!***

!!****f* m_dmfet_driver/dmfet
!! NAME
!! dmfet
!!
!! FUNCTION
!! Primary routine for conducting DMFET calculation.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dmfet(acell,codvsn,dtfil,dtset,mpi_enreg,pawang,pawrad,pawtab,psps,rprim,xred,wvl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet'
!End of the abilint section

 type(dataset_type),intent(inout) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(inout) :: psps
 type(wvl_data),intent(inout) :: wvl !useless
 type(datafiles_type),intent(inout) :: dtfil

 character(len=6),intent(in) :: codvsn
 real(dp),intent(in) :: acell(3),xred(3,dtset%natom)
 real(dp),intent(inout) :: rprim(3,3)

 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawang_type),intent(inout) :: pawang
 
 type(crystal_t) :: crystal
 type(pawfgr_type) :: pawfgr

! type(hdr_type) :: hdr
 type(hdr_type) :: hdr_wfk
! type(ebands_t) :: bstruct 
 type(dmfet_type) :: dmfet_args

 character(len=500) :: message
 integer :: ibtot,isppol,ikibz,ib
 integer :: timrev,option,ierr,comm
 integer :: bantot,mcg,mgfftf,nfftf,my_nspinor,comm_psp,psp_gencond
 integer,parameter :: response=0
 logical :: remove_inv
 real(dp) :: ecut_eff,ecutdg_eff
 real(dp) :: gsqcut_eff,gsqcutc_eff,ucvol
 real(dp) :: ecore

 !arrays
 integer :: npwtot(dtset%nkpt)
 integer,allocatable :: npwarr(:),kg(:,:)
 integer :: ngfft(18),ngfftf(18)
 real(dp) :: gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),allocatable :: cg(:,:)

 character(len=fnlen) :: wfk_fname
 real(dp),pointer :: energies_p(:,:,:)
 real(dp),allocatable :: eigen(:),occ(:)
 real(dp) :: e_fermie

 type(wfk_t) :: wfk
 integer :: band_block(2)
 integer :: iomode, funt
 integer,parameter :: master=0,formeig0=0

!================================================
!The following shouldn't change during DMFET calc
!================================================

 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

 call setup1(acell,bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& dtset%natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,rprim,rprimd,ucvol,psps%usepaw)


!Set up the basis sphere of planewaves
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg, &
& dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
& dtset%mpw,npwarr,npwtot,dtset%nsppol)
 call bandfft_kpt_init1(bandfft_kpt,dtset%istwfk,kg,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol)

 if ( dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   if (psps%useylm==1) then
     option=0
     if (dtset%prtstm==0.and.dtset%iscf>0.and.dtset%positron/=1) option=1 !compute gradients of YLM
     if (dtset%berryopt==4 .and. dtset%optstress /= 0 .and. psps%usepaw==1) option = 1 ! compute gradients of YLM
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
&     psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
   end if
 else
   ABI_ALLOCATE(ylm,(0,0))
   ABI_ALLOCATE(ylmgr,(0,0,0))
 end if

!Open and read pseudopotential files
 comm_psp=mpi_enreg%comm_cell
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,&
& pawrad,pawtab,psps,rprimd,comm_mpi=comm_psp)
! if(psp_gencond==1) write(std_out,*) "psp has been recomputed!"

!========================================
!read wavefunction from previous DFT calc
!========================================
 comm = xmpi_world
 wfk_fname = dtfil%fnamewffk
 call wfk_read_eigenvalues(wfk_fname,energies_p,hdr_wfk,comm)
!check compatability
 call hdr_vs_dtset(hdr_wfk,dtset)

!Crystalline structure
 remove_inv=.false.
 timrev = 2; if (any(dtset%kptopt == [3, 4])) timrev= 1
 call crystal_from_hdr(crystal,hdr_wfk,timrev,remove_inv)
 call crystal_print(crystal)


!FIXME
!gamma point only currently
 ABI_ALLOCATE(eigen,(dtset%mband))
 ABI_ALLOCATE(occ,(dtset%mband))
 eigen(:)=zero; occ(:)=zero

 ibtot=0
 do isppol=1,dtset%nsppol
   do ikibz=1,dtset%nkpt
     do ib=1,hdr_wfk%nband(ikibz+dtset%nkpt*(isppol-1))
       ibtot=ibtot+1
       occ(ibtot)=hdr_wfk%occ(ibtot)
       eigen(ibtot)=energies_p(ib,ikibz,isppol)
     end do
   end do
 end do
 ABI_FREE(energies_p)

 e_fermie = hdr_wfk%fermie

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 mcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol

 ABI_STAT_ALLOCATE(cg,(2,mcg), ierr)
 if(ierr/=0) then
   write(message,"(A)") 'out of memory in cg'
   MSG_ERROR(message)
 endif

 iomode = iomode_from_fname(wfk_fname)
 funt = get_unit()
 do isppol=1,dtset%nsppol
   do ikibz=1,dtset%nkpt
     band_block = [1,hdr_wfk%nband(ikibz+dtset%nkpt*(isppol-1))]
     call wfk_open_read(wfk,wfk_fname,formeig0,iomode,funt,xmpi_comm_self)
     call wfk_read_band_block(wfk,band_block,ikibz,isppol,xmpio_single,cg_k=cg) 
   end do
 end do
!========================================
!finished reading
!========================================

 call dmfet_init(dmfet_args,acell,crystal,dtfil,dtset,psps,mpi_enreg,&
& kg,nfftf,pawfgr,pawtab,pawrad,pawang,npwarr,ylm,ylmgr,mcg,cg,eigen,occ,e_fermie,ecore,wvl)
 call dmfet_run(dmfet_args,rprim,codvsn)
 call destroy_dmfet(dmfet_args)

 call crystal_free(crystal)

 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(cg)

end subroutine dmfet

end module m_dmfet_driver
