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
 use libxc_functionals
 use m_xcdata,           only : get_xclevel,get_auxc_ixc
 use m_mpinfo,           only : distrb2_hf

 use m_crystal,          only : crystal_t,crystal_init,crystal_print
 use defs_wvltypes,      only : wvl_data
 use m_plowannier
 use m_hdr
 use m_ebands
 use m_cgtools

 use m_fft,              only : fourwf,fftpac
 use m_fftcore,          only : sphereboundary
 use m_cgprj,            only : ctocprj
 use m_kg,               only : getph
 use m_pawrhoij,         only : pawrhoij_type,pawrhoij_free,pawrhoij_free_unpacked,symrhoij
 use m_paw_occupancies,  only : initrhoij
 use m_paw_mkrho,        only : denfgr
 use m_paw_mkvij

 use m_pawfgr,           only : pawfgr_type
 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,          only : pawcprj_type, pawcprj_getdim, pawcprj_alloc,pawcprj_free
 use m_pawrad,           only : pawrad_type
 use m_dtset,            only : dtset_copy,dtset_chkneu
 use m_results_gs,       only : results_gs_type,init_results_gs
 use m_ioarr,            only : fftdatar_write
 use m_iowf,             only : outwf

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
   type(pawang_type), pointer :: pawang => null()
   type(pawfgr_type), pointer :: pawfgr => null()

   real(dp), pointer :: acell(:)=>null()
   real(dp), pointer :: ylm(:,:) => null()
   real(dp), pointer :: ylmgr(:,:,:) => null()

   real(dp), pointer :: eigen(:) => null()
   real(dp), pointer :: occ(:) => null()
   real(dp) :: e_fermie,ecore

   complex(dpc),allocatable :: can2sub(:,:)
   real(dp), allocatable ::  occ_wan(:)
   real(dp), allocatable :: sub_occ(:)
   integer :: n_canonical, dim_all, dim_imp,dim_sub, n_frozen

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
& kg,nfftf,pawfgr,pawtab,pawrad,pawang,npwarr,ylm,ylmgr,mcg,cg,eigen,occ,e_fermie,ecore,wvl)


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
 type(pawfgr_type), intent(inout),target :: pawfgr
 type(pawtab_type), intent(inout),target :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrad_type), intent(inout),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type), intent(inout),target :: pawang
 real(dp), intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)

 real(dp), intent(in),target :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

 real(dp), intent(in) :: e_fermie,ecore


 DBG_ENTER("COLL")

 if(dtset%nsubsys.ne.2) MSG_ERROR('Only support nsubsys==2 for now!')

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

 this%pawfgr=>pawfgr
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

 this%pawfgr=>null()
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

 integer :: nband
 integer :: ngfft(18)
 integer :: usecprj,mband_cprj,mcprj,my_nspinor
 integer :: iatom,idir,iorder_cprj,ctocprj_choice,ncpgr
 character(len=500) :: message


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
 !this%dim_all = wan%size_wan
 this%n_canonical = wan%bandf_wan-wan%bandi_wan+1

 call compute_coeff_plowannier(this%crystal,cprj,dimcprj,this%dtset,this%eigen,this%e_fermie,&
&     this%mpi_enreg,this%occ,wan,this%pawtab,this%psps,usecprj,this%dtfil%unpaw,this%pawrad,this%dtfil)


 nband = this%n_canonical
 ABI_ALLOCATE(this%can2sub,(nband,nband))
 ABI_ALLOCATE(this%sub_occ,(nband))
 
! call get_can2sub(wan,1,this%occ,nband,this%can2sub,this%sub_occ,this%dim_imp,this%dim_sub,this%dim_all)

 call dmfet_wan2sub(this,wan)

 write(message,'(3(a,i0,a))') ' dim_imp = ', this%dim_imp, ch10,&
& ' dim_sub = ', this%dim_sub, ch10,&
& ' dim_all = ', this%dim_all, ch10
 call wrtout(std_out,message,'COLL')

! call dmfet_wan2sub(this,wan)

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

 call get_can2sub(wan,onedm_wan,1,this%occ,nbands,this%can2sub,this%sub_occ,this%dim_imp,this%dim_sub,this%dim_all)

! ABI_ALLOCATE(this%occ_wan,(wan%size_wan))
! ABI_ALLOCATE(loc2sub,(wan%size_wan,wan%size_wan))
! call build_subspace(onedm_wan,loc2sub,this%occ_wan,6,wan%size_wan,dim_sub,n_full)
! ABI_ALLOCATE(this%can2sub,(wan%size_wan,wan%size_wan))
! call canonical_to_sub(wan,loc2sub,this%can2sub,wan%size_wan,1,1)


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

! ABI_DEALLOCATE(loc2sub)
 ABI_DEALLOCATE(onedm_wan)

! this%dim_sub = dim_sub
! this%n_frozen = n_full


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

 ABI_ALLOCATE(opersub,(nbands,nbands))
 call compute_oper_ks2sub(operks(1,:,:,1),opersub,this%can2sub,nbands,nbands)

 ABI_ALLOCATE(eig,(nbands))
 ABI_ALLOCATE(rwork,(3*nbands-2))
 lwork = 65*nbands ! Value to optimize speed of the diagonalization
 ABI_ALLOCATE(zwork,(lwork))
 call zheev('v','u',nbands,opersub,nbands,eig,zwork,lwork,rwork,info)
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
subroutine dmfet_core(this,rprim,codvsn)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_core'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 real(dp),intent(inout) :: rprim(3,3)
 character(len=6),intent(in) :: codvsn

 type(dataset_type),allocatable :: sub_dtsets(:)
 type(coeff2_type),allocatable :: dens_sub(:)
 type(oep_type) :: oep_args

 type(results_gs_type) :: res_tot
 type(gstate_sub_input_var) :: scf_inp

 type(hdr_type) :: hdr
 type(ebands_t) :: bstruct

 integer :: nsubsys
 integer :: opt_algorithm = 0
 integer :: dim_sub,i,j
 integer :: me,master

!arrays
 real(dp),allocatable :: dens_tot(:,:),emb_pot(:,:),occ_dummy(:)
 complex(dpc),allocatable::mo_coeff(:,:)

!local MPI
! type(MPI_type) :: l_mpi_enreg


 dim_sub = this%dim_sub

!modify dtset for subspace scf calculation
 this%dtset%mband = this%dim_all
 this%dtset%nband(:) = this%dim_all 
!end modify dtset

! call init_local_mpi_enreg(l_mpi_enreg,this%dtset,this%dtset%mband,this%dtset%nband)

 bstruct = ebands_from_dtset(this%dtset, this%npwarr)
 call hdr_init(bstruct,codvsn,this%dtset,hdr,this%pawtab,0,this%psps,this%wvl%descr,&
& comm_atom=this%mpi_enreg%comm_atom,mpi_atmtab=this%mpi_enreg%my_atmtab)
 call ebands_free(bstruct)

 call print_can2sub(this,this%dtset,this%dtfil,hdr,this%mpi_enreg)

 ABI_ALLOCATE(dens_tot,(dim_sub,dim_sub))

 !total scf calc in subspace
 ABI_ALLOCATE(occ_dummy,(this%dim_all))
 ABI_ALLOCATE(mo_coeff,(this%dim_all,this%dim_all))
 call init_results_gs(this%dtset%natom,this%dtset%nsppol,res_tot)
 call gstate_sub(codvsn,this%acell,this%dtset,this%psps,rprim,res_tot,this%mpi_enreg,this%dtfil,this%wvl,&
& this%cg,this%pawtab,this%pawrad,this%pawang,this%crystal%xred,&
& dens_tot,this%can2sub,this%n_canonical,dim_sub,this%dim_all,this%sub_occ(1:this%dim_all),&
& occ_dummy,mo_coeff,&
& hdr_in=hdr) 

 nsubsys = this%dtset%nsubsys
 ABI_DATATYPE_ALLOCATE(sub_dtsets,(nsubsys))
 call build_subsys(this%dtset,sub_dtsets,nsubsys)

 ABI_ALLOCATE(emb_pot,(dim_sub,dim_sub))
 emb_pot = zero

 ABI_DATATYPE_ALLOCATE(dens_sub,(nsubsys))
 do i=1,nsubsys
   ABI_ALLOCATE(dens_sub(i)%value, (dim_sub,dim_sub))
 enddo

 call gstate_sub_input_var_init(scf_inp,codvsn,this%acell,rprim,this%crystal,&
& this%crystal%xred,this%dtset%natom,this%mcg,this%psps,this%mpi_enreg,&
& this%dtfil,this%wvl,this%cg,this%pawtab,this%pawrad,this%pawang,&
& this%can2sub,this%n_canonical,dim_sub,this%dim_all,this%sub_occ(1:this%dim_all))

 call oep_init(oep_args,scf_inp,dens_tot,dens_sub,emb_pot,opt_algorithm,sub_dtsets,nsubsys)
 call oep_run_split(oep_args,this%dtset%vemb_opt_w_tol,this%dtset%vemb_opt_cycle)


! call print_vemb(this,this%dtset,hdr,this%dtfil,this%crystal,this%mpi_enreg,this%pawfgr,this%kg,this%npwarr,&
!& this%cg,this%mcg,oep_args%V_emb,this%can2sub,this%n_canonical,dim_sub,this%crystal%ucvol)

 call post_energy(this,codvsn,sub_dtsets(1),rprim,oep_args%V_emb,this%n_canonical,this%dim_sub)

 call destroy_oep(oep_args)

!clean memory
 do i=1,nsubsys
   ABI_DEALLOCATE(dens_sub(i)%value)
 enddo
 ABI_DATATYPE_DEALLOCATE(dens_sub)

 ABI_DEALLOCATE(emb_pot)

end subroutine dmfet_core


!!****f* m_dmfet/post_energy
!! NAME
!! post_energy
!!
!! FUNCTION
!! Run a post hybrid functional DFT calculation
!!
!!
!! SOURCE

subroutine post_energy(this,codvsn,dtset,rprim,V_emb,dim_can,dim_sub)

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(dataset_type),intent(inout) :: dtset
 character(len=6),intent(in) :: codvsn
 integer,intent(in) :: dim_can,dim_sub
 real(dp),intent(in) :: V_emb(dim_sub,dim_sub)
 real(dp),intent(inout) :: rprim(3,3)


 type(results_gs_type) :: res_emb
 integer :: ixc_old,ixc
 real(dp),allocatable :: dens_emb(:,:),occ(:),occ_dummy(:)
 complex(dpc),allocatable :: mo_coeff(:,:)

 ABI_ALLOCATE(dens_emb,(dim_sub,dim_sub))
 ABI_ALLOCATE(occ,(dim_sub))
 ABI_ALLOCATE(occ_dummy,(dim_sub))
 ABI_ALLOCATE(mo_coeff,(dim_sub,dim_sub))

 occ =zero; dens_emb=zero

 ixc_old = dtset%ixc
 ixc = -428 !HSE06 
 

 if (ixc_old<0) call libxc_functionals_end()

 dtset%ixc = ixc
!Initialize xclevel and usefock
 call get_xclevel(ixc,dtset%xclevel,dtset%usefock)
 if(dtset%auxc_ixc==0)then
   call get_auxc_ixc(dtset%auxc_ixc,ixc)
 end if
!Now take care of the parameters for hybrid functionals
 if(dtset%usefock==1)then

   if(ixc ==40 .or. ixc ==41 .or. ixc ==42)then
     dtset%hyb_mixing_sr=zero
     dtset%hyb_range_dft=zero ; dtset%hyb_range_fock=zero
     if(ixc==40)dtset%hyb_mixing=one
     if(ixc==41)dtset%hyb_mixing=quarter
     if(ixc==42)dtset%hyb_mixing=third
   else if(ixc==-427)then   ! Special case of HSE03
     dtset%hyb_mixing=zero  ; dtset%hyb_mixing_sr=quarter
     dtset%hyb_range_dft=0.15_dp*two**third  ; dtset%hyb_range_fock=0.15_dp*sqrt(half)
   else if (ixc<0) then
     call libxc_functionals_init(ixc,dtset%nspden)
     call libxc_functionals_get_hybridparams(hyb_mixing=dtset%hyb_mixing,hyb_mixing_sr=dtset%hyb_mixing_sr,&
&     hyb_range=dtset%hyb_range_dft)
     call libxc_functionals_end()
     dtset%hyb_range_fock=dtset%hyb_range_dft
   end if
 endif


 if(dtset%paral_kgb==0)then
   if (dtset%usefock==1) then
       if (dtset%nphf>1) this%mpi_enreg%paral_hf=1
       this%mpi_enreg%nproc_hf = dtset%nphf
       if (dtset%npkpt/=1) then
         this%mpi_enreg%nproc_kpt = dtset%npkpt
       else
         this%mpi_enreg%nproc_kpt = this%mpi_enreg%nproc_cell/this%mpi_enreg%nproc_hf
       end if
   endif
 endif

 if(dtset%paral_kgb>=0) then
     if (dtset%usefock==1) then
       ABI_ALLOCATE(this%mpi_enreg%distrb_hf,(dtset%nkpthf,dtset%nbandhf,1))
!      The dimension of distrb_hf are given by %nkpthf and %nbandhf.
!      We assume that there will be no dependence in spinpol for all the
!      occupied states.
       this%mpi_enreg%distrb_hf=0
     end if

     if(xmpi_paral==1 .and. dtset%usewvl == 0) then
         if (dtset%usefock==1) then
           call distrb2_hf(dtset%nbandhf,dtset%nkpthf,this%mpi_enreg%nproc_cell,dtset%nsppol,this%mpi_enreg)
         end if
     endif
 endif

 if (ixc<0) call libxc_functionals_init(ixc,dtset%nspden)

 write(std_out,*) "Post hybrid functional DFT calculation:"
 write(std_out,*) "usefock = ",dtset%usefock

 call init_results_gs(dtset%natom,dtset%nsppol,res_emb)
 call gstate_sub(codvsn,this%acell,dtset,this%psps,rprim,res_emb,this%mpi_enreg,this%dtfil,this%wvl,&
& this%cg,this%pawtab,this%pawrad,this%pawang,this%crystal%xred,&
& dens_emb,this%can2sub,dim_can,dim_sub,dim_sub,occ,occ_dummy,mo_coeff,&
& emb_pot=V_emb)

 write(std_out,*) "energy of high-level subsystem imp: ",res_emb%etotal

 if (ixc<0) call libxc_functionals_end()
 dtset%ixc = ixc_old
 if (ixc_old<0) call libxc_functionals_init(ixc_old,dtset%nspden)
 

 ABI_DEALLOCATE(dens_emb)
 ABI_DEALLOCATE(occ)

end subroutine post_energy


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
subroutine dmfet_run(this,rprim,codvsn)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_run'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 real(dp),intent(inout)::rprim(3,3)
 character(len=6),intent(in) :: codvsn
 character(len=500) :: message


 write(message,'(2a)') ch10,&
&   '================================================================================='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message,'(2a)') ch10,&
&   '==                           Welcome to sDMFET module                         =='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')


 call dmfet_subspac(this)
 call dmfet_core(this,rprim,codvsn)

end subroutine dmfet_run


!!****f* m_dmfet/build_subsys
!! NAME
!! build_subsys
!!
!! FUNCTION
!! Prepare dtset for subsystems
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
!! dmfet_core
!!
!! CHILDREN
!!
!! SOURCE
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
 character(len=500) :: message

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
  
   if(i==1)then 
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
   else if(i==2)then
    do j=1,dtset%natom
      do k=1,dtset%subsys_natom(1)
        if(j==ia(k)) then
          sub_dtsets(i)%typat(j) = dtset%typat(j) + dtset%ntypat/2
          exit
        endif
      enddo
    enddo
   endif

   call dtset_chkneu(dtset%charge,sub_dtsets(i),dtset%occopt)

   write(message,'(2a,i0,a,f10.6)') ch10,"No. of electrons in subsysetem ",i,": ",sub_dtsets(i)%nelect
   call wrtout(std_out,message,'COLL')

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


!!****f* m_dmfet/print_vemb
!! NAME
!! print_vemb
!!
!! FUNCTION
!! Print embedding potential in real space
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
!! Parallelization
!!
!! PARENTS
!! dmfet_core
!!
!! CHILDREN
!!
!! SOURCE
subroutine print_vemb(this,dtset,hdr,dtfil,crystal,mpi_enreg,pawfgr,kg,npwarr,cg,mcg,vemb,can2sub,norb,nsub,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_vemb'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: crystal
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawfgr_type), intent(in) :: pawfgr

 integer,intent(in) :: norb,nsub,mcg
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)

 real(dp), intent(in) :: ucvol
 real(dp), intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: vemb(nsub,nsub)
 complex(dpc),intent(in) :: can2sub(norb,nsub)

 integer,parameter :: cplex=1,option=0
 integer :: n1,n2,n3,n4,n5,n6,nfft
 integer :: istwf_k,npw_k,nband_k,my_nspinor,tim_fourwf
 integer :: i,j,k, isppol,ikpt,iband, icg,ikg,ipwsp,icwave

 integer :: ngfftf(18),mgfft,ctocprj_choice,iatom,idir,iorder_cprj
 integer :: mcprj,mband_cprj
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 integer,allocatable :: dimcprj_srt(:)

 type(pawcprj_type),pointer :: cprj(:,:)
 type(pawcprj_type),allocatable, target :: cprj_local(:,:)
 type(pawrhoij_type),allocatable :: pawvij(:)
 real(dp),allocatable :: vemb_r(:,:),vemb_r_paw(:,:),vr(:,:,:),vi(:,:,:)
 real(dp),allocatable :: vemb_r_one(:,:),vemb_r_t_one(:,:),nhat_dummy(:,:)
 real(dp),allocatable :: vemb_can_real(:,:),vemb_can_img(:,:)
 real(dp),allocatable :: cwavef(:,:,:),wfraug(:,:,:,:,:)
 real(dp),allocatable :: tmpr(:,:), tmpi(:,:)

 real(dp),allocatable :: ph1d(:,:)

 real(dp) :: dummy(0,0),denpot_dummy(0,0,0)
 
 complex(dpc),allocatable :: tmp(:,:), vemb_cpl(:,:)
 complex(dpc) :: vemb_can(norb,norb)  !vemb in canonical basis, complex or real?

 mgfft=pawfgr%mgfft
 ngfftf(:)=pawfgr%ngfft(:)

!transform vemb to canonical basis
 ABI_ALLOCATE(vemb_cpl,(nsub,nsub))
 do i=1,nsub
   do j=1,nsub
     vemb_cpl(i,j) = cmplx(vemb(i,j),0.0_dp,kind=dp)
   enddo
 enddo

 ABI_ALLOCATE(tmp,(norb,nsub))
 call zgemm('N','N',norb,nsub,nsub,cone,can2sub,norb,vemb_cpl,nsub,czero,tmp,norb)
 call zgemm('N','C',norb,norb,nsub,cone,tmp,norb,can2sub,norb,czero,vemb_can,norb)
 ABI_DEALLOCATE(tmp)
 ABI_DEALLOCATE(vemb_cpl)

! write(std_out,*) " max abs imaginary vemb_can element:", maxval(abs(aimag(vemb_can)))
!end transformation

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 n1 = ngfftf(1) ; n2 = ngfftf(2) ; n3 = ngfftf(3)
 n4 = ngfftf(4) ; n5 = ngfftf(5) ; n6 = ngfftf(6)
 ABI_ALLOCATE(cwavef,(2,dtset%mpw,my_nspinor))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6,norb))

 if(mpi_enreg%paralbd==0) tim_fourwf=3
 if(mpi_enreg%paralbd==1) tim_fourwf=6

 ABI_ALLOCATE(vemb_can_real,(norb,norb))
 ABI_ALLOCATE(vemb_can_img,(norb,norb))
 vemb_can_real = real(vemb_can,kind=dp)
 vemb_can_img = aimag(vemb_can)

 ABI_ALLOCATE(vemb_r,(n1*n2*n3,dtset%nsppol))

 icg = 0
 do isppol=1,dtset%nsppol
   ABI_ALLOCATE(vr,(n4,n5,n6))
   vr = zero
!   ABI_ALLOCATE(vi,(n4,n5,n6))
!   vi = zero

   ikg=0
   do ikpt=1,dtset%nkpt

     nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     npw_k=npwarr(ikpt)
     istwf_k = dtset%istwfk(ikpt)

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg) 
     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     do iband=1,norb
       ipwsp=(iband-1)*npw_k*my_nspinor + icg
       cwavef(:,:,isppol)=cg(:,1+ipwsp:ipwsp+npw_k)

!       write(std_out,*) "debug: <cg|cg>=", iband, cg_zdotc(npw_k,cwavef(:,:,isppol),cwavef(:,:,isppol))

     !cg->ur
       call fourwf(cplex,denpot_dummy,cwavef(:,:,isppol),dummy,wfraug(:,:,:,:,iband),gbound,gbound,istwf_k,&
&        kg_k,kg_k,mgfft,mpi_enreg,1,ngfftf,npw_k,1,n4,n5,n6,option,&
&        dtset%paral_kgb,tim_fourwf,one,one)

!       write(std_out,*) "debug: <ur|ur>=",iband, cg_zdotc(n4*n5*n6,wfraug(:,:,:,:,iband),wfraug(:,:,:,:,iband)),&
!&        cg_zdotc(n4*n5*n6,wfraug(:,:,:,:,iband),wfraug(:,:,:,:,iband))/(n1*n2*n3*1.0d0)
     enddo


     ABI_ALLOCATE(tmpr,(norb,n4*n5*n6))
     ABI_ALLOCATE(tmpi,(norb,n4*n5*n6))
     call dgemm('N','T',norb,n4*n5*n6,norb,one,vemb_can_real,norb,wfraug(1,:,:,:,:),n4*n5*n6,zero,tmpr,norb)
     call dgemm('N','T',norb,n4*n5*n6,norb,one,vemb_can_img,norb,wfraug(2,:,:,:,:),n4*n5*n6,one,tmpr,norb)
     call dgemm('N','T',norb,n4*n5*n6,norb,one,vemb_can_img,norb,wfraug(1,:,:,:,:),n4*n5*n6,zero,tmpi,norb)
     call dgemm('N','T',norb,n4*n5*n6,norb,-one,vemb_can_real,norb,wfraug(2,:,:,:,:),n4*n5*n6,one,tmpi,norb)


     do i=1,n6
       do j=1,n5
         do k=1,n4
           vr(k,j,i) = vr(k,j,i) + dot_product(wfraug(1,k,j,i,:), tmpr(:,k+(j-1)*n4+(i-1)*n5*n4) )
           vr(k,j,i) = vr(k,j,i) - dot_product(wfraug(2,k,j,i,:), tmpi(:,k+(j-1)*n4+(i-1)*n5*n4) )
!           vi(k,j,i) = vi(k,j,i) + dot_product(wfraug(1,k,j,i,:), tmpi(:,k+(j-1)*n4+(i-1)*n5*n4) )
!           vi(k,j,i) = vi(k,j,i) + dot_product(wfraug(2,k,j,i,:), tmpr(:,k+(j-1)*n4+(i-1)*n5*n4) )
         enddo
       enddo
     enddo

     ABI_DEALLOCATE(tmpr)
     ABI_DEALLOCATE(tmpi)

     icg=icg+npw_k*my_nspinor*nband_k
     ikg=ikg+npw_k
   enddo !ikpt

   call fftpac(isppol,mpi_enreg,dtset%nsppol,n1,n2,n3,n4,n5,n6,ngfftf,vemb_r,vr,1)
   vemb_r = vemb_r/ucvol
!   write(std_out,*) "max(vemb_r_aug)=",maxval(vr)
!   write(std_out,*) "max(vemb_r)=",maxval(vemb_r)

!   write(std_out,*) " max abs imaginary vemb_r element:", maxval(abs(vi))

   ABI_DEALLOCATE(vr)
!   ABI_DEALLOCATE(vi)
 enddo !isppol


!PAW part
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 call getph(this%crystal%atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,this%crystal%xred)

!compute cprj
 ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
 call pawcprj_getdim(dimcprj_srt,dtset%natom,this%crystal%nattyp,dtset%ntypat,dtset%typat,this%pawtab,'O')
 mband_cprj=dtset%mband
 if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
 mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
 ABI_DATATYPE_ALLOCATE(cprj_local,(dtset%natom,mcprj))
 ctocprj_choice = 1
 call pawcprj_alloc(cprj_local,0,dimcprj_srt)
 iatom=0 ; iorder_cprj=0 !0: sorted; 1:un-sorted
 idir = 0
 call ctocprj(this%crystal%atindx,this%cg,ctocprj_choice,cprj_local,this%crystal%gmet,this%crystal%gprimd,&
&   iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&   mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,this%psps%mpsang,&
&   dtset%mpw,dtset%natom,this%crystal%nattyp,dtset%nband,dtset%natom,dtset%ngfft,&
&   dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,dtset%paral_kgb,ph1d,this%psps,this%crystal%rmet,dtset%typat,&
&   ucvol,dtfil%unpaw,this%crystal%xred,this%ylm,this%ylmgr)

 ABI_DEALLOCATE(ph1d)

 cprj=> cprj_local

 ABI_DATATYPE_ALLOCATE(pawvij,(dtset%natom))
 call initrhoij(dtset%pawcpxocc,dtset%lexexch,&
&   dtset%lpawu,dtset%natom,dtset%natom,dtset%nspden,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,pawvij,dtset%pawspnorb,this%pawtab,dtset%spinat,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 call pawmkvij(this%crystal%atindx,this%crystal%atindx1,vemb_can,norb,cprj,dimcprj_srt,&
&   dtset%istwfk,dtset%mband,mband_cprj,mcprj,dtset%mkmem,mpi_enreg,&
&   dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,dtset%paral_kgb,&
&   pawvij,dtfil%unpaw)

 ABI_DEALLOCATE(dimcprj_srt)

!symmetrize potential and store it in packed form
 call symrhoij(pawvij,pawvij,1,this%crystal%gprimd,this%crystal%indsym,0,dtset%natom,dtset%nsym,dtset%ntypat,&
&   1,this%pawang,dtset%pawprtvol,this%pawtab,this%crystal%rprimd,dtset%symafm,this%crystal%symrec,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 call pawrhoij_free_unpacked(pawvij)

!full potential on fine grid
 ABI_ALLOCATE(vemb_r_paw,(pawfgr%nfft,dtset%nspden))
 ABI_ALLOCATE(vemb_r_one,(pawfgr%nfft,dtset%nspden))
 ABI_ALLOCATE(vemb_r_t_one,(pawfgr%nfft,dtset%nspden))
 ABI_ALLOCATE(nhat_dummy,(pawfgr%nfft,dtset%nspden))
 vemb_r_paw = zero
 vemb_r_one = zero
 vemb_r_t_one = zero
 nhat_dummy = zero !no use
 call denfgr(this%crystal%atindx1,this%crystal%gmet,mpi_enreg%comm_fft,dtset%natom,dtset%natom,&
&   this%crystal%nattyp,ngfftf,nhat_dummy,dtset%nspinor,dtset%nsppol,dtset%nspden,&
&   dtset%ntypat,pawfgr,this%pawrad,pawvij,this%pawtab,dtset%prtvol,vemb_r,vemb_r_paw,vemb_r_one,&
&   vemb_r_t_one,this%crystal%rprimd,dtset%typat,ucvol,this%crystal%xred,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!clean up
 call pawcprj_free(cprj_local)
 ABI_DATATYPE_DEALLOCATE(cprj_local)

 call pawrhoij_free(pawvij)
 ABI_DATATYPE_DEALLOCATE(pawvij)

! ABI_DEALLOCATE(vemb_r)
 ABI_DEALLOCATE(vemb_r_one)
 ABI_DEALLOCATE(vemb_r_t_one) 
 ABI_DEALLOCATE(nhat_dummy)

!**********************
!*  write to disk     *
!**********************
 nfft = n1*n2*n3
 call fftdatar_write("vemb",dtfil%fnameabo_app_vemb,dtset%iomode,hdr,&
&     crystal,ngfftf,cplex,nfft,dtset%nspden,vemb_r,mpi_enreg)

!clean up
 ABI_DEALLOCATE(vemb_r_paw)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(vemb_can_real)
 ABI_DEALLOCATE(vemb_can_img)

end subroutine print_vemb


subroutine print_can2sub(this,dtset,dtfil,hdr,mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_can2sub'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(datafiles_type),intent(in) :: dtfil
 type(MPI_type),intent(in) :: mpi_enreg

 integer :: dim_can,dim_all,mcg
 integer :: nband_sub(1)

 real(dp),allocatable::cg_new(:,:)
 real(dp),allocatable::tmp_real(:,:),tmp_img(:,:)
 real(dp),allocatable::resid(:),eigen(:)

 dim_can = this%n_canonical
 dim_all =this%dim_all

 ABI_ALLOCATE(tmp_real,(dim_can,dim_all))
 ABI_ALLOCATE(tmp_img,(dim_can,dim_all))
 tmp_real=real(this%can2sub)
 tmp_img=aimag(this%can2sub)

 mcg = dtset%mpw*dim_all
 ABI_ALLOCATE(cg_new,(2,mcg))
 call dgemm('N','N',dtset%mpw,dim_all,dim_can,one,this%cg(1,:),dtset%mpw,tmp_real,dim_can,zero,cg_new(1,:),dtset%mpw)
 call dgemm('N','N',dtset%mpw,dim_all,dim_can,-one,this%cg(2,:),dtset%mpw,tmp_img,dim_can,one,cg_new(1,:),dtset%mpw)
 call dgemm('N','N',dtset%mpw,dim_all,dim_can,one,this%cg(1,:),dtset%mpw,tmp_img,dim_can,zero,cg_new(2,:),dtset%mpw)
 call dgemm('N','N',dtset%mpw,dim_all,dim_can,one,this%cg(2,:),dtset%mpw,tmp_real,dim_can,one,cg_new(2,:),dtset%mpw)

 nband_sub = dim_all
 hdr%nband = nband_sub
 hdr%mband = dim_all
 ABI_ALLOCATE(resid,(dim_all*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(eigen,(dim_all*dtset%nkpt*dtset%nsppol))
 eigen(:)=zero ; resid(:)=zero
 call outwf(cg_new,dtset,this%psps,eigen,dtfil%fnameabo_suborb,hdr,this%kg,dtset%kptns,&
&   dim_all,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,&
&   nband_sub,dtset%nkpt,this%npwarr,dtset%nsppol,&
&   this%sub_occ(1:dim_all),resid,0,dtfil%unwff2,this%wvl%wfs,this%wvl%descr)

end subroutine print_can2sub


end module m_dmfet
