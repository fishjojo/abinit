!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dmfet
!! NAME
!!  m_dmfet
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


module m_dmfet

 use defs_basis
 use m_abicore
 use m_errors
 use defs_datatypes
 use defs_abitypes

 use m_crystal,          only : crystal_t,crystal_init,crystal_print
 use defs_wvltypes,      only : wvl_data
 use m_plowannier

 use m_cgprj,            only : ctocprj
 use m_kg,               only : getph
 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,          only : pawcprj_type, pawcprj_getdim, pawcprj_alloc,pawcprj_free
 use m_pawrad,           only : pawrad_type
 use m_dtset,            only : dtset_copy,dtset_chkneu
 use m_results_gs,       only : results_gs_type,init_results_gs

 use m_dmfet_oep
 use m_gstate_sub
 use m_subscf

 implicit none

 private

 public :: dmfet_init
 public :: dmfet_run
 public :: destroy_dmfet


 type, public :: dmfet_type

   integer,pointer :: nfftf => null()
   integer,pointer :: mcg => null()
   real(dp), pointer :: cg(:,:) => null()

   type(datafiles_type),pointer :: dtfil => null()
   type(dataset_type),pointer :: dtset => null()
   type(pseudopotential_type),pointer :: psps => null()
   type(crystal_t),pointer :: crystal => null()
   type(MPI_type),pointer :: mpi_enreg => null()

   integer, pointer :: kg(:,:) => null()
   integer, pointer :: npwarr(:) => null()
   type(pawtab_type), pointer :: pawtab(:) => null()
   type(pawrad_type), pointer :: pawrad(:) => null()
   type(pawang_type),pointer :: pawang => null()

   real(dp), pointer :: acell(:)=>null()
   real(dp), pointer :: ylm(:,:) => null()
   real(dp), pointer :: ylmgr(:,:,:) => null()

   real(dp), pointer :: eigen(:) => null()
   real(dp), pointer :: occ(:) => null()
   real(dp) :: e_fermie,ecore

   complex(dpc),allocatable :: can2sub(:,:)
   real(dp), allocatable ::  occ_wan(:)
   integer :: n_canonical, dim_all, dim_sub, n_frozen

   type(wvl_data),pointer :: wvl => null()

 end type dmfet_type


contains
!!****f* m_dmfet/dmfet_init
!! NAME
!! dmfet_init
!!
!! FUNCTION
!! Initialize a dmfet object
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
!! m_dmfet_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine dmfet_init(this,acell,crystal,dtfil,dtset,psps,mpi_enreg,&
& kg,nfftf,pawtab,pawrad,pawang,npwarr,ylm,ylmgr,mcg,cg,eigen,occ,e_fermie,ecore,wvl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_init'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(datafiles_type),intent(in),target :: dtfil
 type(dataset_type),intent(in),target :: dtset
 type(pseudopotential_type),intent(in),target :: psps
 type(crystal_t),intent(in),target :: crystal
 type(MPI_type),intent(in),target :: mpi_enreg
 type(wvl_data),intent(in),target :: wvl !useless
 integer, intent(in),target :: npwarr(dtset%nkpt)
 real(dp),intent(in),target :: acell(3)


 integer,intent(in),target :: nfftf,mcg
 real(dp), intent(in),target :: cg(2,mcg)
 integer, intent(in),target :: kg(3,dtset%mpw*dtset%mkmem)
 type(pawtab_type), intent(inout),target :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrad_type), intent(inout),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type), intent(inout),target :: pawang
 real(dp), intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)

 real(dp), intent(in),target :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

 real(dp), intent(in) :: e_fermie,ecore

 logical :: DEBUG=.FALSE.

 DBG_ENTER("COLL")

 this%dtfil=>dtfil
 this%dtset=>dtset
 this%psps=>psps
 this%acell=>acell
 this%crystal=>crystal
 this%mpi_enreg=>mpi_enreg
 this%nfftf=>nfftf

 this%kg=>kg
 this%npwarr=>npwarr
 this%mcg=>mcg
 this%cg=>cg

 this%pawtab=>pawtab
 this%pawrad=>pawrad
 this%pawang=>pawang

 this%ylm=>ylm
 this%ylmgr=>ylmgr

 this%wvl=>wvl

 this%eigen=>eigen
 this%occ=>occ
 this%e_fermie = e_fermie
 this%ecore = ecore

 this%n_canonical = 0
 this%dim_all = 0
 this%dim_sub = 0
 this%n_frozen = 0

 DBG_EXIT("COLL")

end subroutine dmfet_init


!!****f* m_dmfet/destroy_dmfet
!! NAME
!! destroy_dmfet
!!
!! FUNCTION
!! Destroy a dmfet object
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!! m_dmfet_driver
!!
!! SOURCE
subroutine destroy_dmfet(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_dmfet'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this

 this%dtfil=>null()
 this%dtset=>null()
 this%psps=>null()
 this%acell=>null()
 this%crystal=>null()
 this%mpi_enreg=>null()
 this%nfftf=>null()

 this%kg=>null()
 this%npwarr=>null()
 this%mcg=>null()
 this%cg=>null()

 this%pawtab=>null()
 this%pawrad=>null()
 this%pawang=>null()

 this%ylm=>null()
 this%ylmgr=>null()

 this%wvl=>null()

 this%eigen=>null()
 this%occ=>null()

 this%e_fermie = -1.d5
 this%ecore = -1.d5

 this%n_canonical = 0
 this%dim_all = 0
 this%dim_sub = 0
 this%n_frozen = 0

 if(allocated(this%can2sub)) ABI_DEALLOCATE(this%can2sub)
 if(allocated(this%occ_wan)) ABI_DEALLOCATE(this%occ_wan)

end subroutine destroy_dmfet

!!****f* m_dmfet/dmfet_subspac
!! NAME
!! dmfet_subspac
!!
!! FUNCTION
!! Construct the subspace for a dmfet calculation
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
subroutine dmfet_subspac(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_subspac'
!End of the abilint section

 type(dmfet_type), intent(inout) :: this

 type(plowannier_type) :: wan
 type(pawcprj_type),allocatable, target :: cprj_local(:,:)
 type(pawcprj_type),pointer :: cprj(:,:)
 integer,allocatable :: dimcprj(:)
 real(dp),allocatable :: ph1d(:,:)


 integer :: ngfft(18)
 integer :: usecprj,mband_cprj,mcprj,my_nspinor
 integer :: iatom,idir,iorder_cprj,ctocprj_choice,ncpgr


 ABI_ALLOCATE(dimcprj,(this%dtset%natom))
 call pawcprj_getdim(dimcprj,this%dtset%natom,this%crystal%nattyp,this%dtset%ntypat,this%dtset%typat,this%pawtab,'R')

 ngfft(:)=this%dtset%ngfft(:)
 ABI_ALLOCATE(ph1d,(2,3*(2*this%dtset%mgfft+1)*this%dtset%natom))
 call getph(this%crystal%atindx,this%dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,this%crystal%xred)

 usecprj=1
 my_nspinor=max(1,this%dtset%nspinor/this%mpi_enreg%nproc_spinor)
 mband_cprj=this%dtset%mband/this%mpi_enreg%nproc_band
 mcprj=my_nspinor*mband_cprj*this%dtset%mkmem*this%dtset%nsppol
 ABI_DATATYPE_ALLOCATE(cprj_local,(this%dtset%natom,mcprj))
 ncpgr = 0 ; ctocprj_choice = 1
 call pawcprj_alloc(cprj_local,ncpgr,dimcprj)
 cprj=> cprj_local
 iatom=0 ; iorder_cprj=1 ! cprj are not ordered
 idir = 0
 call ctocprj(this%crystal%atindx,this%cg,ctocprj_choice,cprj_local,this%crystal%gmet,this%crystal%gprimd,&
&   iatom,idir,iorder_cprj,this%dtset%istwfk,this%kg,this%dtset%kptns,&
&   this%mcg,mcprj,this%dtset%mgfft,this%dtset%mkmem,this%mpi_enreg,this%psps%mpsang,&
&   this%dtset%mpw,this%dtset%natom,this%crystal%nattyp,this%dtset%nband,this%dtset%natom,ngfft,&
&   this%dtset%nkpt,this%dtset%nloalg,this%npwarr,this%dtset%nspinor,this%dtset%nsppol,&
&   this%dtset%ntypat,this%dtset%paral_kgb,ph1d,this%psps,this%crystal%rmet,this%dtset%typat,&
&   this%crystal%ucvol,this%dtfil%unpaw,this%crystal%xred,this%ylm,this%ylmgr)


 call init_plowannier(this%dtset,wan)
 this%dim_all = wan%size_wan
 this%n_canonical = wan%bandf_wan-wan%bandi_wan+1

 call compute_coeff_plowannier(this%crystal,cprj,dimcprj,this%dtset,this%eigen,this%e_fermie,&
&     this%mpi_enreg,this%occ,wan,this%pawtab,this%psps,usecprj,this%dtfil%unpaw,this%pawrad,this%dtfil)

 call dmfet_wan2sub(this,wan)

 call destroy_plowannier(wan)

 ABI_DEALLOCATE(ph1d)
 call pawcprj_free(cprj_local)
 ABI_DATATYPE_DEALLOCATE(cprj_local)
 ABI_DEALLOCATE(dimcprj)

end subroutine dmfet_subspac


!!****f* m_dmfet/dmfet_wan2sub
!! NAME
!! dmfet_wan2sub
!!
!! FUNCTION
!! Construct subspace orbitals using Wannier orbitals
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
subroutine dmfet_wan2sub(this,wan)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_wan2sub'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(plowannier_type),intent(inout) :: wan

 type(operwan_type), allocatable :: operwan(:,:,:)
 complex(dpc), allocatable :: operks(:,:,:,:)

 real(dp), allocatable :: onedm_wan(:,:,:,:), loc2sub(:,:)
 !complex(dpc), allocatable :: uuT(:,:,:), uTu(:,:,:)
 complex(dpc), allocatable :: opersub(:,:)
 integer:: dim_sub, n_full

 integer :: nbands,isppol,iband1,ibandc,ikpt
 integer :: ii,jj

 real(dp), allocatable :: eig(:), rwork(:)
 complex(dpc), allocatable :: zwork(:)
 integer :: lwork,info



 nbands = this%n_canonical

 !Inialize an empty Wannier operator
 ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
 call initialize_operwan(wan,operwan)

 !Creation of the KS occupation operator
 ABI_ALLOCATE(operks,(wan%nkpt,nbands,nbands,wan%nsppol))
 operks = czero
 do isppol = 1,wan%nsppol
   do iband1 = 1,nbands
     ibandc = iband1 + wan%bandi_wan - 1
     do ikpt = 1,wan%nkpt
       operks(ikpt,iband1,iband1,isppol) = this%occ(((ikpt-1)*this%dtset%mband+ibandc+(isppol-1)*wan%nkpt*this%dtset%mband))
     end do
   end do
 end do

 !compute the occupation in wannier basis
 do ikpt = 1,wan%nkpt
   call compute_oper_ks2wan(wan,operks,operwan,ikpt)
 end do
 ABI_DEALLOCATE(operks)

 ABI_ALLOCATE(onedm_wan,(wan%size_wan,wan%size_wan,wan%nkpt,wan%nsppol))
 call oper_to_matrix(wan,operwan,onedm_wan)

 call destroy_operwan(wan,operwan)
 ABI_DATATYPE_DEALLOCATE(operwan)
!   write(std_out,*)"debug:1pdm_wan"
!   do ii=1,wan%size_wan
!      write(std_out,*) onedm_wan(ii,:,1,1) 
!   enddo

 ABI_ALLOCATE(this%occ_wan,(wan%size_wan))
 ABI_ALLOCATE(loc2sub,(wan%size_wan,wan%size_wan))
 call build_subspace(onedm_wan,loc2sub,this%occ_wan,6,wan%size_wan,dim_sub,n_full)
 ABI_ALLOCATE(this%can2sub,(wan%size_wan,wan%size_wan))
 call canonical_to_sub(wan,loc2sub,this%can2sub,wan%size_wan,1,1)


!   ABI_ALLOCATE(uuT,(wan%size_wan,wan%size_wan,wan%nkpt))
!   ABI_ALLOCATE(uTu,(wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nkpt))
!   call test_unitary(wan,uuT,uTu,1)
!   write(std_out,*)"debug:uuT"
!   do ii=1,wan%size_wan
!      write(std_out,*) uuT(ii,:,1)
!   enddo

!   write(std_out,*)"debug:uTu"
!   do ii=1,wan%bandf_wan-wan%bandi_wan+1
!      write(std_out,*) uTu(ii,:,1)
!   enddo

 ABI_DEALLOCATE(loc2sub)
 ABI_DEALLOCATE(onedm_wan)

 this%dim_sub = dim_sub
 this%n_frozen = n_full


 ! KS eigenvalue matrix
 ABI_ALLOCATE(operks,(wan%nkpt,nbands,nbands,wan%nsppol))
 operks = czero
 do isppol = 1,wan%nsppol
   do iband1 = 1,nbands
     ibandc = iband1 + wan%bandi_wan - 1
     do ikpt = 1,wan%nkpt
       operks(ikpt,iband1,iband1,isppol) = this%eigen(((ikpt-1)*this%dtset%mband+ibandc+(isppol-1)*wan%nkpt*this%dtset%mband))
     end do
   end do
 end do

 ABI_ALLOCATE(opersub,(wan%size_wan,wan%size_wan))
 call compute_oper_ks2sub(operks(1,:,:,1),opersub,this%can2sub,nbands,wan%size_wan)

 ABI_ALLOCATE(eig,(wan%size_wan))
 ABI_ALLOCATE(rwork,(3*wan%size_wan-2))
 lwork = 65*wan%size_wan ! Value to optimize speed of the diagonalization
 ABI_ALLOCATE(zwork,(lwork))
 call zheev('v','u',wan%size_wan,opersub,wan%size_wan,eig,zwork,lwork,rwork,info)
 write(ab_out,*) "debug: eigen_sub"
 write(ab_out,*) eig(:)

 ABI_DEALLOCATE(operks)
 ABI_DEALLOCATE(eig)
 ABI_DEALLOCATE(rwork)
 ABI_DEALLOCATE(zwork)
 ABI_DEALLOCATE(opersub)

end subroutine dmfet_wan2sub



!!****f* m_dmfet/dmfet_core
!! NAME
!! dmfet_core
!!
!! FUNCTION
!! Main function for performing dmfet calculations using subspace orbitals
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
subroutine dmfet_core(this,rprim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_core'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 real(dp),intent(inout) :: rprim(3,3)

 type(dataset_type),allocatable :: sub_dtsets(:)
 type(coeff2_type),allocatable :: dens_sub(:)
 type(oep_type) :: oep_args

 type(results_gs_type) :: res_tot
 type(gstate_sub_input_var) :: scf_inp

 integer :: nsubsys
 integer :: opt_algorithm = 0
 integer :: dim_sub,i,j

!arrays
 real(dp),allocatable :: dens_tot(:,:),emb_pot(:,:)

 dim_sub = this%dim_all !use all of sub orbitals for now
 ABI_ALLOCATE(dens_tot,(dim_sub,dim_sub))

 !total scf calc in subspace
 call init_results_gs(this%dtset%natom,this%dtset%nsppol,res_tot)
 call gstate_sub(this%acell,this%dtset,this%psps,rprim,res_tot,this%mpi_enreg,this%dtfil,this%wvl,&
& this%cg,this%pawtab,this%pawrad,this%pawang,this%crystal%xred,&
& dens_tot,this%can2sub,this%dim_all,dim_sub) 


 nsubsys = this%dtset%nsubsys
 ABI_DATATYPE_ALLOCATE(sub_dtsets,(nsubsys))
 call build_subsys(this%dtset,sub_dtsets,nsubsys)

 ABI_ALLOCATE(emb_pot,(dim_sub,dim_sub))
 emb_pot = zero

 ABI_DATATYPE_ALLOCATE(dens_sub,(nsubsys))
 do i=1,nsubsys
   ABI_ALLOCATE(dens_sub(i)%value, (dim_sub,dim_sub))
 enddo

 call gstate_sub_input_var_init(scf_inp,this%acell,rprim,this%crystal%xred,this%dtset%natom,this%mcg,this%psps,this%mpi_enreg,&
& this%dtfil,this%wvl,this%cg,this%pawtab,this%pawrad,this%pawang,this%can2sub,this%dim_all,dim_sub)

 call oep_init(oep_args,scf_inp,dens_tot,dens_sub,emb_pot,opt_algorithm,sub_dtsets,nsubsys)
 call oep_run(oep_args,this%dtset%vemb_opt_w_tol,this%dtset%vemb_opt_cycle)
 call destroy_oep(oep_args)


!clean memory
 do i=1,nsubsys
   ABI_DEALLOCATE(dens_sub(i)%value)
 enddo
 ABI_DATATYPE_DEALLOCATE(dens_sub)

end subroutine dmfet_core



!!****f* m_dmfet/dmfet_run
!! NAME
!! dmfet_run
!!
!! FUNCTION
!! Run a dmfet calculation
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
!! m_dmfet_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine dmfet_run(this,rprim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_run'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 real(dp),intent(inout)::rprim(3,3)

 character(len=500) :: message


 write(message,'(2a,i3)') ch10,&
&   '================================================================================='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message,'(2a,i3)') ch10,&
&   '==                           Welcome to sDMFET module                         =='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')


 call dmfet_subspac(this)
 call dmfet_core(this,rprim)

end subroutine dmfet_run



subroutine build_subsys(dtset,sub_dtsets,nsubsys)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'build_subsys'
!End of the abilint section

 implicit none

 integer,intent(in) :: nsubsys
 type(dataset_type),intent(in) :: dtset
 type(dataset_type),intent(inout) :: sub_dtsets(nsubsys)
! type(crystal_t),intent(in) :: crystal
! type(crystal_t),intent(inout) :: sub_crystals(nsubsys)
 !type(coeff2_type),intent(inout) :: sub_xreds(nsubsys)

 integer :: i, j, k, ioff, sub_natom
 integer,allocatable :: ia(:)

 ABI_ALLOCATE(ia,(dtset%natom))

 ioff = 0
 do i=1,nsubsys
   call dtset_copy(sub_dtsets(i), dtset)
   sub_natom = dtset%subsys_natom(i)
   !sub_dtsets(i)%natom = sub_natom

   !ABI_DEALLOCATE(sub_dtsets(i)%typat)
   !ABI_ALLOCATE(sub_dtsets(i)%typat,(sub_natom))

   !ABI_DEALLOCATE(sub_dtsets(i)%xred_orig)
   !ABI_ALLOCATE(sub_dtsets(i)%xred_orig,(3,sub_natom,1))

   !ABI_ALLOCATE(sub_xreds(i)%value,(3,sub_natom))
   
   do j=1,sub_natom
     ia(j) = dtset%subsys_iatom(j+ioff)
     !sub_dtsets(i)%typat(j) = dtset%typat(ia)
     !sub_dtsets(i)%xred_orig(:,j,1) = dtset%xred_orig(:,ia,1)
     !sub_xreds(i)%value(:,j) = crystal%xred(:,ia)
   enddo

   do j=1,dtset%natom
      do k=1,sub_natom
        if(ia(k).eq.j) goto 200
      enddo
      sub_dtsets(i)%typat(j) = dtset%typat(j) + dtset%ntypat/2
200 enddo

   call dtset_chkneu(dtset%charge,sub_dtsets(i),dtset%occopt)
   write(std_out,*) "No. of electrons in subsysetem ",i,": ",sub_dtsets(i)%nelect 
   ioff = ioff + sub_natom
 enddo


! do i=1,nsubsys
!   call crystal_init(crystal%amu,sub_crystals(i),crystal%space_group,dtset%natom,&
!&   crystal%npsp,crystal%ntypat,crystal%nsym,crystal%rprimd,sub_dtsets(i)%typat,crystal%xred,&
!&   crystal%zion,crystal%znucl,crystal%timrev,crystal%use_antiferro,.false.,crystal%title,&
!&   crystal%symrel,crystal%tnons,crystal%symafm)
! enddo

 ABI_DEALLOCATE(ia)

end subroutine build_subsys


end module m_dmfet
