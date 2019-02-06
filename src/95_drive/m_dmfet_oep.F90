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

 use m_crystal,          only : crystal_t
 use m_gstate_sub

 implicit none

 private

 public :: oep_init
 public :: oep_run
 public :: destroy_oep


 type, public :: oep_type

   real(dp),pointer :: P_ref(:,:)=>null()
   real(dp),pointer :: V_emb(:,:)=>null()
   complex(dpc),pointer :: can2sub(:,:)=>null()
   type(coeff2_type),pointer :: P_sub(:)=>null()
   type(dataset_type),pointer :: sub_dtsets(:)=>null()

   integer,pointer :: dim_can=>null()
   integer,pointer :: dimP=>null()
   integer,pointer :: opt_algorithm=>null()
   integer,pointer :: nsubsys=>null()

 end type oep_type


contains

subroutine oep_init(this,P_ref,P_sub,V_emb,can2sub,dim_can,dimP,opt_algorithm,sub_dtsets,nsubsys)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'oep_init'
!End of the abilint section

 implicit none

 type(oep_type),intent(inout) :: this

 integer,intent(in),target :: dim_can,dimP,opt_algorithm,nsubsys
 real(dp),intent(in),target :: P_ref(dimP,dimP)
 real(dp),intent(inout),target :: V_emb(dimP,dimP)
 complex(dpc),intent(in),target :: can2sub(dim_can,dimP)
 type(coeff2_type),intent(inout),target :: P_sub(nsubsys)
 type(dataset_type),intent(inout),target :: sub_dtsets(nsubsys)

 this%dim_can=>dim_can
 this%dimP=>dimP
 this%nsubsys=>nsubsys
 this%opt_algorithm=>opt_algorithm


 this%P_ref=>P_ref
 this%V_emb=>V_emb
 this%P_sub=>P_sub
 this%can2sub=>can2sub

 this%sub_dtsets=>sub_dtsets

end subroutine oep_init




subroutine oep_run(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'oep_run'
!End of the abilint section

 implicit none

 include 'nlopt.f'

 type(oep_type),intent(inout) :: this

 integer*8 opt  !pointer to the opt object
 integer algorithm, n, ires
 integer i,j,k,dimP

 real(dp) :: minf, tol
 real(dp),allocatable :: x(:)

 opt = 0
 algorithm = NLOPT_LD_LBFGS !temporarily hard coded
 dimP = this%dimP
 n = dimP*(dimP+1)/2

 ABI_ALLOCATE(x,(n))
 call mat2vec(x,this%V_emb,dimP)

 call nlo_create(opt, algorithm, n)

 tol = 1.d-5
 call nlo_set_ftol_abs(ires, opt, tol)

 call nlo_set_max_objective(ires, opt, cost_wuyang, this)

 call nlo_optimize(ires, opt, x, minf)

 call nlo_destroy(opt)

end subroutine oep_run


subroutine mat2vec(vec,mat,n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mat2vec'
!End of the abilint section

 implicit none

 integer,intent(in) :: n
 real(dp),intent(in) :: mat(n,n)
 real(dp),intent(out) :: vec(n*(n+1)/2)
 integer :: i,j,k

 k=1
 do j=1,n
   do i=j,n
     vec(k) = mat(i,j) !lower triangle
     k=k+1
   enddo
 enddo

end subroutine mat2vec


subroutine vec2mat(vec,mat,n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vec2mat'
!End of the abilint section

 implicit none

 integer,intent(in) :: n
 real(dp),intent(out) :: mat(n,n)
 real(dp),intent(in) :: vec(n*(n+1)/2)
 integer :: i,j,k

 k=1
 do j=1,n
   do i=j,n
     mat(i,j) = vec(k) !lower triangle
     mat(j,i) = mat(i,j)
     k=k+1
   enddo
 enddo

end subroutine vec2mat


function trace_dot(A,B,n) result(trace)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'trace_dot'
!End of the abilint section

 implicit none

 integer,intent(in) :: n
 real(dp),intent(in) :: A(n,n), B(n,n)
 real(dp) :: trace
 integer :: i,j

 trace = 0.0_dp
 do i=1,n
   do j=1,n
     trace = trace + A(i,j)*B(j,i)
   enddo
 enddo 

end function trace_dot


subroutine cost_wuyang(f, n, x, grad, need_gradient, this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cost_wuyang'
!End of the abilint section

 implicit none
     
 integer n,need_gradient
 integer i, dim_can, dimP
 real(dp) :: f, x(n), grad(n), e_imp, e_bath
 real(dp),allocatable :: diffP(:,:), e_sub(:)
 type(oep_type) :: this

 dim_can = this%dim_can
 dimP = this%dimP
 call vec2mat(x,this%V_emb,dimP)

 !subsystem scf
 ABI_ALLOCATE(e_sub,(this%nsubsys))
! do i=1,this%nsubsys
!   call gstate_sub(this%sub_crystals(i),this%sub_dtsets(i),this%psps,this%mpi_enreg,this%dtfil,&
!&   this%npwarr,this%nfftf,this%mcg,this%cg,this%ylm,this%ylmgr,this%kg,&
!&   this%pawtab,this%pawrad,this%pawang,this%pawfgr,this%wvl,&
!&   this%P_sub(i),this%can2sub,dim_can,dimP,this%ecore,this%occ,e_sub(i),this%V_emb)
! enddo


 !objective function
 f = 0.0_dp
 do i=1,this%nsubsys
   f = f + e_sub(i)
 enddo 
 f = f - trace_dot(this%P_ref,this%V_emb,dimP)

 ABI_DEALLOCATE(e_sub)

 !gradient
 if (need_gradient.ne.0) then
   !compute gradient
   ABI_ALLOCATE(diffP,(dimP,dimP))
   diffP = zero
   do i=1,this%nsubsys
     diffP = diffP + this%P_sub(i)%value
   enddo
   diffP = diffP - this%P_ref
   call mat2vec(grad,diffP,dimP)
   ABI_DEALLOCATE(diffP)
 endif

 
end subroutine cost_wuyang



subroutine destroy_oep(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_oep'
!End of the abilint section

 implicit none
 
 type(oep_type),intent(inout) :: this


 this%dimP=>null()
 this%nsubsys=>null()
 this%opt_algorithm=>null()

 this%P_ref=>null()
 this%V_emb=>null()
 this%P_sub=>null()

 this%sub_dtsets=>null()


end subroutine destroy_oep

end module m_dmfet_oep

