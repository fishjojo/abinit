!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_oep
!! NAME
!!  m_oep
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


module m_dmfet_oep

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore

 use m_gstate_sub

 implicit none

 private

 public :: oep_init
 public :: oep_run
 public :: destroy_oep


 type, public :: oep_type

   real(dp),pointer :: P_ref(:,:)=>null()
   real(dp),pointer :: P_imp(:,:)=>null()
   real(dp),pointer :: P_bath(:,:)=>null()
   real(dp),pointer :: V_emb(:,:)=>null()

   integer,pointer :: dimP=>null()
   integer,pointer :: opt_algorithm=>null()

 end type oep_type


contains

subroutine oep_init(this,P_ref,V_emb,dimP,opt_algorithm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'oep_init'
!End of the abilint section

 implicit none

 type(oep_type),intent(inout) :: this

 integer,intent(in),target :: dimP,opt_algorithm
 real(dp),intent(in),target :: P_ref(dimP,dimP)
 real(dp),intent(inout),target :: V_emb(dimP,dimP)

 this%dimP=>dimP

 this%opt_algorithm=>opt_algorithm

 this%P_ref=>P_ref
 this%V_emb=>V_emb


end subroutine oep_init




subroutine oep_run(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'oep_run'
!End of the abilint section

 implicit none

 type(oep_type),intent(inout) :: this


end subroutine oep_run




subroutine destroy_oep(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_oep'
!End of the abilint section

 implicit none
 
 type(oep_type),intent(inout) :: this

 this%dimP=>null()

 this%opt_algorithm=>null()

 this%P_ref=>null()
 this%V_emb=>null()


end subroutine destroy_oep

end module m_dmfet_oep

