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
 use m_errors

 use m_crystal,          only : crystal_t
 use m_gstate_sub
 use m_results_gs,       only : results_gs_type,init_results_gs,destroy_results_gs

 implicit none

 private

 public :: oep_init
 public :: oep_run
 public :: oep_run_split
 public :: destroy_oep


 type, public :: oep_type

   type(gstate_sub_input_var),pointer :: scf_inp=>null()

   integer,pointer :: opt_algorithm=>null()
   integer,pointer :: nsubsys=>null()

   real(dp),pointer :: P_ref(:,:)=>null()
   real(dp),pointer :: V_emb(:,:)=>null()

   type(coeff2_type),pointer :: P_sub(:)=>null()
   type(dataset_type),pointer :: sub_dtsets(:)=>null()

   real(dp),allocatable :: occ_imp(:),occ_bath(:)
   complex(dpc),allocatable :: mo_coeff_imp(:,:), mo_coeff_bath(:,:)
   complex(dpc),allocatable :: fock_imp(:,:),fock_bath(:,:)


   !debug
!   type(hdr_type),pointer :: hdr=>null()
!   type(crystal_t),pointer :: crystal_tot => null()

 end type oep_type


contains

subroutine oep_init(this,scf_inp,P_ref,P_sub,V_emb,opt_algorithm,sub_dtsets,nsubsys)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'oep_init'
!End of the abilint section

 implicit none

 type(oep_type),intent(inout) :: this
 type(gstate_sub_input_var),intent(in),target :: scf_inp

 integer,intent(in),target :: opt_algorithm,nsubsys
 real(dp),intent(in),target :: P_ref(scf_inp%dim_sub,scf_inp%dim_sub)
 real(dp),intent(inout),target :: V_emb(scf_inp%dim_sub,scf_inp%dim_sub)
 type(coeff2_type),intent(inout),target :: P_sub(nsubsys)
 type(dataset_type),intent(inout),target :: sub_dtsets(nsubsys)

 !debug
! type(hdr_type),intent(inout),target :: hdr
! type(crystal_t),intent(in),target ::crystal_tot

! this%hdr => hdr
! this%crystal_tot=>crystal_tot

 this%scf_inp=>scf_inp
 this%nsubsys=>nsubsys
 this%opt_algorithm=>opt_algorithm

 this%P_ref=>P_ref
 this%V_emb=>V_emb
 this%P_sub=>P_sub

 this%sub_dtsets=>sub_dtsets

end subroutine oep_init



subroutine oep_run_split(this,wtol,max_cycle)

 implicit none

 type(oep_type),intent(inout) :: this
 integer,intent(in) :: wtol,max_cycle

 integer :: dim_all_loc, dim_sub, dim_all,dim_can
 integer :: i,iter,nstep

 real(dp):: gtol,dist_imp,dist_bath
 character(len=500)::message
 type(results_gs_type),allocatable :: results(:)
 real(dp),allocatable :: P_imp_old(:,:),P_bath_old(:,:)

 dim_can = this%scf_inp%dim_can
 dim_sub = this%scf_inp%dim_sub
 dim_all = this%scf_inp%dim_all

 ABI_ALLOCATE(this%occ_imp,(dim_sub))
 ABI_ALLOCATE(this%occ_bath,(dim_all))
 ABI_ALLOCATE(this%mo_coeff_imp,(dim_sub,dim_sub))
 ABI_ALLOCATE(this%mo_coeff_bath,(dim_all,dim_all))
 ABI_ALLOCATE(this%fock_imp,(dim_sub,dim_sub))
 ABI_ALLOCATE(this%fock_bath,(dim_sub,dim_sub))

 ABI_DATATYPE_ALLOCATE(results,(this%nsubsys))

 do i=1,this%nsubsys
   call init_results_gs(this%sub_dtsets(i)%natom,this%sub_dtsets(i)%nsppol,results(i))
   if(i==1)then
     !embedded region should not have frozen orbitals
     dim_all_loc = dim_sub
   else
     dim_all_loc = dim_all
   endif
   this%sub_dtsets(i)%mband = dim_all_loc
   this%sub_dtsets(i)%nband = dim_all_loc
   if(i==1)then
    call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(i),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(i),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(i)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_all_loc,this%scf_inp%sub_occ,&
&   this%occ_imp,this%mo_coeff_imp,&
&   emb_pot=this%V_emb)
   else
    call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(i),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(i),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(i)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_all_loc,this%scf_inp%sub_occ,&
&   this%occ_bath,this%mo_coeff_bath,&
&   emb_pot=this%V_emb)
   endif

   !e_sub(i) = results(i)%etotal
   !write(message,'(2a,i0,a,f16.10)') ch10,"energy of subsystem ",i,": ",e_sub(i)
   !call wrtout(std_out,message,'COLL')
 enddo

 nstep = this%sub_dtsets(1)%nstep
 this%sub_dtsets(1)%nstep = 1
 this%sub_dtsets(2)%nstep = 1

 ABI_ALLOCATE(P_imp_old,(dim_sub,dim_sub))
 ABI_ALLOCATE(P_bath_old,(dim_sub,dim_sub))
 gtol = 1.0d-5

 do iter=1,max_cycle

   write(message,'(2a,i0)') ch10,"OEP ITER ", iter
   call wrtout(std_out,message,'COLL')

!  compute Fock matrix
   !call init_results_gs(this%sub_dtsets(1)%natom,this%sub_dtsets(1)%nsppol,results(1))
   call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(1),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(1),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(1)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_sub,this%scf_inp%sub_occ,&
&   this%occ_imp,this%mo_coeff_imp,&
&   emb_pot=this%V_emb,fock_mat=this%fock_imp)

   !call destroy_results_gs(results(1))

   !call init_results_gs(this%sub_dtsets(2)%natom,this%sub_dtsets(2)%nsppol,results(2))
   call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(2),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(2),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(2)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_all,this%scf_inp%sub_occ,&
&   this%occ_bath,this%mo_coeff_bath,&
&   emb_pot=this%V_emb,fock_mat=this%fock_bath)

   !call destroy_results_gs(results(2))

   P_imp_old = this%P_sub(1)%value
   P_bath_old = this%P_sub(2)%value

!  minimize W with fixed Fock
   call opt_vemb_fix_fock(this,wtol,max_cycle)

!  check convergence
   dist_imp = maxval(abs(this%P_sub(1)%value - P_imp_old))
   dist_bath = maxval(abs(this%P_sub(2)%value - P_bath_old))
   write(std_out,*) 'dist_imp=',dist_imp,'dist_bath=',dist_bath
   if(dist_imp<gtol.and.dist_bath<gtol) then
     exit
   endif

 enddo

 do i=1,this%nsubsys
   call destroy_results_gs(results(i))
 enddo
 ABI_DATATYPE_DEALLOCATE(results)

 this%sub_dtsets(1)%nstep = nstep
 this%sub_dtsets(2)%nstep = nstep
 call verify_scf(this)

end subroutine oep_run_split


subroutine verify_scf(this)

 implicit none

 type(oep_type),intent(inout) :: this
 integer :: dim_can,dim_sub,dim_all,i
 type(results_gs_type),allocatable :: results(:)
 real(dp),allocatable :: diffP(:,:) 

 dim_can = this%scf_inp%dim_can
 dim_sub = this%scf_inp%dim_sub
 dim_all = this%scf_inp%dim_all

 ABI_DATATYPE_ALLOCATE(results,(this%nsubsys))

 do i=1,this%nsubsys
   call init_results_gs(this%sub_dtsets(i)%natom,this%sub_dtsets(i)%nsppol,results(i))
 enddo

 call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(1),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(1),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(1)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_sub,this%scf_inp%sub_occ,&
&   this%occ_imp,this%mo_coeff_imp,&
&   emb_pot=this%V_emb,fock_mat=this%fock_imp,prtden=1,&
&   crystal_tot=this%scf_inp%crystal_tot)
 write(std_out,*) "energy of subsystem imp: ",results(1)%etotal

 call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(2),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(2),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(2)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_all,this%scf_inp%sub_occ,&
&   this%occ_bath,this%mo_coeff_bath,&
&   emb_pot=this%V_emb,fock_mat=this%fock_bath,prtden=2,&
&   crystal_tot=this%scf_inp%crystal_tot)
 write(std_out,*) "energy of subsystem bath: ",results(2)%etotal

 ABI_ALLOCATE(diffP,(dim_sub,dim_sub))
 diffP = zero
 do i=1,this%nsubsys
   diffP = diffP + this%P_sub(i)%value
 enddo
 diffP = diffP - this%P_ref

 write(std_out,*) 'max element of diffP:', maxval(abs(diffP))
 ABI_DEALLOCATE(diffP)

 do i=1,this%nsubsys
   call destroy_results_gs(results(i))
 enddo
 ABI_DATATYPE_DEALLOCATE(results)

end subroutine verify_scf


subroutine opt_vemb_fix_fock(this,wtol,max_cycle)

 implicit none

 include 'nlopt.f'

 type(oep_type),intent(inout) :: this
 integer,intent(in) :: wtol,max_cycle

 integer*8 opt  !pointer to the opt object
 integer algorithm, n, ires
 integer dim_sub,maxeval
 real(dp)::tol,minf
 real(dp),allocatable :: x(:)

 opt = 0
 algorithm = NLOPT_LD_LBFGS !temporarily hard coded
 !algorithm = NLOPT_LD_MMA
 dim_sub = this%scf_inp%dim_sub
 n = dim_sub*(dim_sub+1)/2

 ABI_ALLOCATE(x,(n))
 call mat2vec(x,this%V_emb,dim_sub)

 call nlo_create(opt, algorithm, n)

 tol = 10.0_dp**(-wtol)
 call nlo_set_ftol_rel(ires, opt, tol)
 call nlo_get_ftol_rel(tol, opt)
 write(std_out,*) 'debug: tol=',tol

 call nlo_set_maxeval(ires, opt, max_cycle)
 call nlo_get_maxeval(maxeval, opt)
 write(std_out,*) 'debug: maxeval',maxeval

 call nlo_set_min_objective(ires, opt, cost_wuyang_fix_fock, this)

 call nlo_optimize(ires, opt, x, minf)

 if(ires<0) MSG_WARNING('nlopt failed')

 call nlo_destroy(opt)

end subroutine opt_vemb_fix_fock




subroutine cost_wuyang_fix_fock(f, n, x, grad, need_gradient, this)

 implicit none

 integer n,need_gradient
 integer i,ii, dim_can, dim_sub, dim_all, dim_all_loc
 real(dp) :: f, x(n), grad(n), e_imp, e_bath
 real(dp),allocatable :: diffP(:,:), e_sub(:)
 type(oep_type) :: this

 type(results_gs_type),allocatable :: results(:)
 real(dp),allocatable :: eig_imp(:),eig_bath(:)
 complex(dpc),allocatable :: ham_imp(:,:),ham_bath(:,:),dens_mat(:,:),tmp(:,:)

 dim_can = this%scf_inp%dim_can
 dim_sub = this%scf_inp%dim_sub
 dim_all = this%scf_inp%dim_all
 call vec2mat(x,this%V_emb,dim_sub)

 !write(std_out,*) 'Vemb:'
 !write(std_out,*) this%V_emb


 ABI_ALLOCATE(ham_imp,(dim_sub,dim_sub))
 ABI_ALLOCATE(eig_imp,(dim_sub))
 ham_imp = this%fock_imp + this%V_emb

 !write(std_out,*) 'ham_imp'
 !write(std_out,*) ham_imp
 call diag_hermit_mat(ham_imp,eig_imp,dim_sub)
 this%mo_coeff_imp = ham_imp

 ABI_ALLOCATE(ham_bath,(dim_sub,dim_sub))
 ABI_ALLOCATE(eig_bath,(dim_sub))
 ham_bath = this%fock_bath + this%V_emb
 call diag_hermit_mat(ham_bath,eig_bath,dim_sub)
 this%mo_coeff_bath(1:dim_sub,1:dim_sub) = ham_bath

 !imp
 !write(std_out,*) 'occ_imp'
 !write(std_out,*) this%occ_imp
 ABI_ALLOCATE(dens_mat,(dim_sub,dim_sub))
 dens_mat=czero
 do ii=1,dim_sub
   dens_mat(ii,ii) = cmplx(this%occ_imp(ii),kind=dp)
 enddo

 ABI_ALLOCATE(tmp,(dim_sub,dim_sub))
 call zgemm('N','C',dim_sub,dim_sub,dim_sub,cone,dens_mat,dim_sub,ham_imp,dim_sub,czero,tmp,dim_sub)
 call zgemm('N','N',dim_sub,dim_sub,dim_sub,cone,ham_imp,dim_sub,tmp,dim_sub,czero,dens_mat,dim_sub)
 ABI_DEALLOCATE(tmp)

 if(maxval(abs(aimag(dens_mat)))>tol8) then
   MSG_WARNING('Density matrix is not real!')
   write(std_out,*) "max abs imaginary element:", maxval(abs(aimag(dens_mat)))
 endif

 e_imp = real(trace_dot_cmplx(this%fock_imp+this%V_emb,dens_mat,dim_sub),kind=dp)
 this%P_sub(1)%value = real(dens_mat,kind=dp)
 ABI_DEALLOCATE(dens_mat)


 !bath
 !write(std_out,*) 'occ_bath'
 !write(std_out,*) this%occ_bath
 ABI_ALLOCATE(dens_mat,(dim_sub,dim_sub))
 dens_mat=czero
 do ii=1,dim_sub
   dens_mat(ii,ii) = cmplx(this%occ_bath(ii),kind=dp)
 enddo

 ABI_ALLOCATE(tmp,(dim_sub,dim_sub))
 call zgemm('N','C',dim_sub,dim_sub,dim_sub,cone,dens_mat,dim_sub,ham_bath,dim_sub,czero,tmp,dim_sub)
 call zgemm('N','N',dim_sub,dim_sub,dim_sub,cone,ham_bath,dim_sub,tmp,dim_sub,czero,dens_mat,dim_sub)
 ABI_DEALLOCATE(tmp)

 if(maxval(abs(aimag(dens_mat)))>tol8) then
   MSG_WARNING('Density matrix is not real!')
   write(std_out,*) "max abs imaginary element:", maxval(abs(aimag(dens_mat)))
 endif

 e_bath = real(trace_dot_cmplx(this%fock_bath+this%V_emb,dens_mat,dim_sub),kind=dp)
 this%P_sub(2)%value = real(dens_mat,kind=dp)
 ABI_DEALLOCATE(dens_mat)

 !objective function
 f = e_imp + e_bath
 f = f - trace_dot(this%P_ref,this%V_emb,dim_sub)

 write(std_out,*) "objective function value:", f

 f = -1.0_dp*f


 !gradient
 if (need_gradient.ne.0) then
   !compute gradient
   ABI_ALLOCATE(diffP,(dim_sub,dim_sub))
   diffP = zero
   do i=1,this%nsubsys
     diffP = diffP + this%P_sub(i)%value
   enddo
   diffP = diffP - this%P_ref
   !write(std_out,*) 'diffP'
   !write(std_out,*) diffP
   call mat2vec(grad,diffP,dim_sub)
   grad = -1.0_dp*grad
   write(std_out,*) "max component of gradient:", maxval(abs(grad))
   ABI_DEALLOCATE(diffP)
 endif

end subroutine cost_wuyang_fix_fock



subroutine diag_hermit_mat(mat,w,dim)

 implicit none
 
 integer,intent(in)::dim
 real(dp),intent(inout):: w(dim)
 complex(dpc),intent(inout) :: mat(dim,dim)

 integer :: lwork,info
 real(dp), allocatable :: rwork(:)
 complex(dpc), allocatable :: zwork(:)

 ABI_ALLOCATE(rwork,(3*dim-2))
 lwork = 65*dim ! Value to optimize speed of the diagonalization
 ABI_ALLOCATE(zwork,(lwork))
 call zheev('v','u',dim,mat,dim,w,zwork,lwork,rwork,info)
 if(info.ne.0) MSG_ERROR('Hamiltonian diagonalization failed!')

 ABI_DEALLOCATE(zwork)
 ABI_DEALLOCATE(rwork)

end subroutine diag_hermit_mat


subroutine oep_run(this,wtol,max_cycle)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'oep_run'
!End of the abilint section

 implicit none

 include 'nlopt.f'

 type(oep_type),intent(inout) :: this
 integer,intent(in) :: wtol,max_cycle

 integer*8 opt  !pointer to the opt object
 integer algorithm, n, ires
 integer i,j,k,dim_sub,maxeval

 real(dp) :: minf, tol
 real(dp),allocatable :: x(:)

 opt = 0
 !algorithm = NLOPT_LD_LBFGS !temporarily hard coded
 algorithm = NLOPT_LD_MMA
 dim_sub = this%scf_inp%dim_sub
 n = dim_sub*(dim_sub+1)/2

 ABI_ALLOCATE(x,(n))
 call mat2vec(x,this%V_emb,dim_sub)

 call nlo_create(opt, algorithm, n)

 tol = 10.0_dp**(-wtol)
 call nlo_set_ftol_rel(ires, opt, tol)
 call nlo_get_ftol_rel(tol, opt)
 write(std_out,*) 'debug: tol=',tol

 call nlo_set_maxeval(ires, opt, max_cycle)
 call nlo_get_maxeval(maxeval, opt)
 write(std_out,*) 'debug: maxeval',maxeval

 call nlo_set_min_objective(ires, opt, cost_wuyang, this)

 call nlo_optimize(ires, opt, x, minf)

 if(ires<0) MSG_WARNING('nlopt failed')

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

function trace_dot_cmplx(A,B,n) result(trace)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'trace_dot'
!End of the abilint section

 implicit none

 integer,intent(in) :: n
 complex(dpc),intent(in) :: A(n,n), B(n,n)
 complex(dpc) :: trace
 integer :: i,j

 trace = czero
 do i=1,n
   do j=1,n
     trace = trace + A(i,j)*B(j,i)
   enddo
 enddo

end function trace_dot_cmplx



subroutine cost_wuyang(f, n, x, grad, need_gradient, this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cost_wuyang'
!End of the abilint section

 implicit none
     
 integer n,need_gradient
 integer i, dim_can, dim_sub, dim_all, dim_all_loc
 real(dp) :: f, x(n), grad(n), e_imp, e_bath
 real(dp),allocatable :: diffP(:,:), e_sub(:)
 type(oep_type) :: this

 type(results_gs_type),allocatable :: results(:)
 real(dp),allocatable :: occ(:)
 complex(dpc),allocatable :: mo_coeff(:,:)

 dim_can = this%scf_inp%dim_can
 dim_sub = this%scf_inp%dim_sub
 dim_all = this%scf_inp%dim_all
 call vec2mat(x,this%V_emb,dim_sub)

! write(std_out,*) "Vemb:"
! do i=1,dim_sub
!   write(std_out,*) this%V_emb(i,:)
! enddo

 !subsystem scf
 ABI_ALLOCATE(e_sub,(this%nsubsys))

 ABI_DATATYPE_ALLOCATE(results,(this%nsubsys))

 do i=1,this%nsubsys
   call init_results_gs(this%sub_dtsets(i)%natom,this%sub_dtsets(i)%nsppol,results(i))
   if(i==1)then
     !embedded region should not have frozen orbitals
     dim_all_loc = dim_sub
   else
     dim_all_loc = dim_all
   endif
   this%sub_dtsets(i)%mband = dim_all_loc
   this%sub_dtsets(i)%nband = dim_all_loc
   ABI_ALLOCATE(occ,(dim_all_loc))
   ABI_ALLOCATE(mo_coeff,(dim_all_loc,dim_all_loc))
   call gstate_sub(this%scf_inp%codvsn,this%scf_inp%acell,this%sub_dtsets(i),this%scf_inp%psps,this%scf_inp%rprim,&
&   results(i),this%scf_inp%mpi_enreg,this%scf_inp%dtfil,this%scf_inp%wvl,&
&   this%scf_inp%cg,this%scf_inp%pawtab,this%scf_inp%pawrad,this%scf_inp%pawang,this%scf_inp%xred,&
&   this%P_sub(i)%value,this%scf_inp%can2sub,dim_can,dim_sub,dim_all_loc,this%scf_inp%sub_occ,&
&   occ,mo_coeff,&
&   emb_pot=this%V_emb)

   e_sub(i) = results(i)%etotal
   write(std_out,*) "energy of subsystem ",i,": ",e_sub(i)

   ABI_DEALLOCATE(occ)
   ABI_DEALLOCATE(mo_coeff)
 enddo


 !objective function
 f = 0.0_dp
 do i=1,this%nsubsys
   f = f + e_sub(i)
 enddo 
 f = f - trace_dot(this%P_ref,this%V_emb,dim_sub)

 write(std_out,*) "objective function value:", f

 f = -1.0_dp*f

 ABI_DEALLOCATE(e_sub)

 !gradient
 if (need_gradient.ne.0) then
   !compute gradient
   ABI_ALLOCATE(diffP,(dim_sub,dim_sub))
   diffP = zero
   do i=1,this%nsubsys
     diffP = diffP + this%P_sub(i)%value
   enddo
   diffP = diffP - this%P_ref
   call mat2vec(grad,diffP,dim_sub)
   grad = -1.0_dp*grad
   write(std_out,*) "max component of gradient:", maxval(abs(grad))
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


 this%nsubsys=>null()
 this%nsubsys=>null()
 this%opt_algorithm=>null()

 this%P_ref=>null()
 this%V_emb=>null()
 this%P_sub=>null()

 this%sub_dtsets=>null()


end subroutine destroy_oep

end module m_dmfet_oep

