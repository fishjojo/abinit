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



subroutine cost_wuyang(f, n, x, grad, need_gradient, f_data)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cost_wuyang'
!End of the abilint section

 implicit none
     
 integer n,need_gradient
 integer i
 real(dp) :: f, x(n), grad(n)
 type(oep_type),intent(inout) :: f_data

 call vec2mat(x,f_data%V_emb,f_data%dimP)


 !imp scf

 !bath scf

 !f=e_imp+e_bath

 if (need_gradient.ne.0) then
   !compute gradient

   !grad = P_imp + P_bath - P_ref
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

 this%opt_algorithm=>null()

 this%P_ref=>null()
 this%V_emb=>null()


end subroutine destroy_oep

end module m_dmfet_oep

