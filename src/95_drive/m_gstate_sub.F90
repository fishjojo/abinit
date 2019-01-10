!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gstate_sub
!! NAME
!!  m_gstate_sub
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


module m_gstate_sub
 use defs_basis
 use m_abicore
 use defs_datatypes
 use defs_abitypes

 use m_crystal,          only : crystal_t
 use defs_wvltypes,      only : wvl_data

 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_pawfgr,           only : pawfgr_type
 use m_pawrad,           only : pawrad_type

 use m_subscf

 implicit none

 private

 public :: gstate_sub

contains

subroutine gstate_sub(crystal,dtset,psps,mpi_enreg,dtfil,&
& npwarr,nfftf,mcg,cg,ylm,ylmgr,kg,&
& pawtab,pawrad,pawang,pawfgr,wvl,&
& dens,can2sub,dim_can,dim_sub,ecore,occ,&
& emb_pot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gstate_sub'
!End of the abilint section

 implicit none

 type(crystal_t),intent(in) :: crystal
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil

 integer, intent(in) :: npwarr(dtset%nkpt)
 integer, intent(in) :: nfftf,mcg
 real(dp), intent(in) :: cg(2,mcg)

 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)

 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrad_type), intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type),intent(inout) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr

 type(wvl_data),intent(inout) :: wvl !useless

 integer, intent(in) :: dim_can, dim_sub
 real(dp), intent(in), optional :: emb_pot(dim_sub,dim_sub)
 complex(dpc), intent(in) :: can2sub(dim_can,dim_sub)

 real(dp), intent(inout) :: ecore
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: dens(dim_sub,dim_sub)

 type(subscf_type) :: subscf_args
!arrays

 if(present(emb_pot)) then
  call subscf_init(subscf_args,dtfil,dtset,psps,crystal,&
&  nfftf,pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&  ylm,ylmgr,kg,cg,mcg,npwarr,ecore,wvl,occ,0,dens,dim_sub,emb_pot)
 else
  call subscf_init(subscf_args,dtfil,dtset,psps,crystal,&
&  nfftf,pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&  ylm,ylmgr,kg,cg,mcg,npwarr,ecore,wvl,occ,0,dens,dim_sub)
 endif

 call subscf_run(subscf_args,can2sub,dim_can,dim_sub)

 call subscf_destroy(subscf_args)

end subroutine gstate_sub






end module m_gstate_sub
