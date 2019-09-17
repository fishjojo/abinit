!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_subscf
!! NAME
!!  m_subscf
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

module m_subscf

 use defs_basis
 use m_abicore
 use m_errors
 use defs_datatypes
 use defs_abitypes

 use m_xg
 use m_cgtools
 use m_efield
 use m_gemm_nonlop
 use m_hdr
 use m_xmpi

 use m_ab7_mixing

 use m_cgwf,             only : mksubham
 use m_spacepar,         only : meanvalue_g
 use m_scfcv_core,       only : etotfor,wf_mixing
 use m_rhotov,           only : rhotov
 use m_common,           only : scprqt 
 use m_newrho,           only : newrho
 use m_occ,              only : newocc
 use m_plowannier,       only : compute_oper_ks2sub
 use defs_wvltypes,      only : wvl_data
 use m_geometry,         only : metric,fixsym
 use m_crystal,          only : crystal_t
 use m_kg,               only : getph,getcut,mkkin,mkkpg
 use m_mkffnl,           only : mkffnl
 use m_pawtab,           only : pawtab_type,pawtab_get_lsize
 use m_paw_an,           only : paw_an_type,paw_an_init,paw_an_free,paw_an_nullify,paw_an_reset_flags
 use m_paw_ij,           only : paw_ij_type,paw_ij_init,paw_ij_free,paw_ij_nullify,paw_ij_reset_flags
 use m_pawfgr,           only : pawfgr_type
 use m_pawrad,           only : pawrad_type
 use m_pawang,           only : pawang_type
 use m_pawfgrtab,        only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_paw_nhat,         only : nhatgrid,pawmknhat
 use m_paw_tools,        only : chkpawovlp
 use m_pawrhoij,         only : pawrhoij_type
 use m_paw_occupancies,  only : initrhoij
 use m_paw_denpot,       only : pawdenpot
 use m_pawdij,           only : pawdij, symdij,pawdijhat
 use m_energies,         only : energies_type, energies_init, energies_copy
 use m_scf_history,      only : scf_history_type, scf_history_init, scf_history_free
 use m_fock,             only : fock_type, fock_init, fock_destroy, fock_ACE_destroy, fock_common_destroy, &
                                fock_BZ_destroy, fock_update_exc, fock_updatecwaveocc, fock_set_ieigen,fock_updateikpt
 use m_fock_getghc,      only : fock2ACE,fock_ACE_getghc

 use m_symtk,            only : symmetrize_xred
 use m_spacepar,         only : setsym
 use m_drivexc,          only : check_kxc
 use m_setvtr,           only : setvtr
 use m_mkrho,            only : initro,mkrho,prtrhomxmn
 use m_hamiltonian,      only : init_hamiltonian,destroy_hamiltonian, &
&                               load_spin_hamiltonian, load_k_hamiltonian,gs_hamiltonian_type
 use m_electronpositron, only : electronpositron_type

 use m_fourier_interpol, only : transgrid
 use m_fft,              only : fftpac,fourdp
 use m_getghc,           only : getghc,multithreaded_getghc

 use m_paw_dmft,         only : paw_dmft_type
 use m_cgprj,            only : ctocprj
 use m_paw_occupancies,  only : pawmkrhoij
 use m_paw_mkrho,        only : pawmkrho
 use m_pawcprj,          only : pawcprj_type, pawcprj_alloc,pawcprj_getdim,pawcprj_free

 use m_results_gs,       only : results_gs_type
 use m_prep_kgb,         only : prep_getghc
 use m_bandfft_kpt
 use m_mpinfo,           only : proc_distrb_cycle

 !debug
 use m_ebands
 use m_ioarr,            only : fftdatar_write

 implicit none

 private 

 public :: subscf_init
 public :: subscf_run
 public :: subscf_destroy
 public :: cgtosub

!!****t* m_subscf/subscf_type
!! NAME
!!  subscf_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE
 type, public :: subscf_type

   type(datafiles_type),pointer :: dtfil => null()
   type(dataset_type),pointer :: dtset => null()
   type(crystal_t),pointer :: crystal => null()
   type(MPI_type),pointer :: mpi_enreg => null()
   type(pseudopotential_type),pointer :: psps => null()
   type(wvl_data),pointer :: wvl => null()

   type(pawtab_type), pointer :: pawtab(:) => null()
   type(pawfgr_type),pointer :: pawfgr => null()
   type(pawang_type),pointer :: pawang => null()
   type(pawrad_type), pointer :: pawrad(:) => null()
   type(pawrhoij_type),pointer :: pawrhoij(:) => null()
   type(pawcprj_type),pointer :: cprj(:,:) => null()

   integer, pointer :: kg(:,:) => null()
   integer, pointer :: npwarr(:) => null()
   real(dp),pointer :: ylm(:,:) => null()
   real(dp),pointer :: ylmgr(:,:,:) => null()

   integer :: ireadwf
   integer, pointer :: nfftf => null()
   integer, pointer :: mcg => null()
   integer,pointer :: mcprj => null()
   integer, pointer :: my_natom => null()
   real(dp),pointer :: ecore => null()
   real(dp),pointer :: cg(:,:) => null()
   real(dp),pointer :: occ(:) => null()

   real(dp),pointer :: rhog(:,:) => null()
   real(dp),pointer :: rhor(:,:) => null()

   real(dp),allocatable :: eig_sub(:)
   complex(dpc),pointer :: subham_sub(:,:) => null()
   complex(dpc),pointer :: fock_mat(:,:) => null()

   real(dp), pointer :: dens_mat_real(:,:,:) => null()
   real(dp), pointer :: emb_pot(:,:,:) => null()

   type(results_gs_type),pointer :: results_gs => null()
   logical :: has_embpot,save_fock_mat

   !useless
   type(electronpositron_type),pointer :: electronpositron => null()
   type(paw_dmft_type),pointer :: paw_dmft => null()
   type(efield_type),pointer :: dtefield => null()
   integer, pointer :: pwind_alloc => null()
   integer, pointer :: pwind(:,:,:) => null()
   real(dp),pointer :: pwnsfac(:,:) => null()
   real(dp),pointer :: taug(:,:)=>null()
   real(dp),pointer :: taur(:,:)=>null()

 end type subscf_type

contains 

!!****f* m_subscf/subscf_init
!! NAME
!! subscf_init
!!
!! FUNCTION
!! Initialize a subscf object
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
subroutine subscf_init(this,dtfil,dtset,psps,results_gs,crystal,nfftf,&
&                      pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&                      ylm,ylmgr,kg,cg,mcg,cprj,mcprj,my_natom,npwarr,ecore,wvl,&
&                      occ,rhog,rhor,pawrhoij,dens_mat_real,dim_sub,dim_all,&
&                      taug,taur,paw_dmft,dtefield,pwind_alloc,pwind,pwnsfac,electronpositron,& !useless
&                      emb_pot,ireadwf,fock_mat,subham_sub) !optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_init'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout):: this
 type(datafiles_type),intent(in),target :: dtfil
 type(dataset_type),intent(in),target :: dtset
 type(crystal_t),intent(in),target :: crystal
 type(MPI_type),intent(in),target :: mpi_enreg
 type(pseudopotential_type),intent(in),target :: psps
 type(pawtab_type), intent(in),target :: pawtab(psps%ntypat*psps%usepaw)
 type(pawfgr_type),intent(in),target :: pawfgr
 type(pawrad_type), intent(in),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type),intent(in),target :: pawang
 type(wvl_data),intent(in),target :: wvl !useless
 type(results_gs_type),intent(inout),target :: results_gs

 integer, intent(in),target :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in),target :: nfftf,mcg,mcprj,my_natom
 integer, intent(in),target :: npwarr(dtset%nkpt)
 integer, intent(in),target,optional :: ireadwf
 integer, intent(in) :: dim_sub,dim_all


 real(dp),intent(in),target :: occ(dim_all*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout),target :: cg(2,mcg),ecore
 real(dp),intent(in),target :: rhor(nfftf,dtset%nspden),rhog(2,nfftf)
 real(dp),intent(in),target :: dens_mat_real(2,dim_sub,dim_sub)
 real(dp),intent(in),target, optional :: emb_pot(2,dim_sub,dim_sub)
 complex(dpc),intent(inout),target,optional :: fock_mat(dim_sub,dim_sub),subham_sub(dim_all,dim_all)

 type(pawrhoij_type), intent(in),target :: pawrhoij(my_natom*psps%usepaw)
 type(pawcprj_type), allocatable,intent(in),target :: cprj(:,:)

!useless
 type(paw_dmft_type),intent(in),target :: paw_dmft
 type(efield_type),intent(in),target :: dtefield
 type(electronpositron_type),target :: electronpositron  
 integer, intent(in),target :: pwind_alloc
 integer, intent(in),target :: pwind(pwind_alloc,2,3)
 real(dp),intent(in),target :: pwnsfac(2,pwind_alloc)
 real(dp),intent(in),target :: taug(2,nfftf*dtset%usekden)
 real(dp),intent(in),target :: taur(nfftf,dtset%nspden*dtset%usekden)

 integer :: ii

! write(std_out,*) "Entering subscf_init"
 !initialize pointers
 this%dtfil=>dtfil
 this%dtset=>dtset
 this%psps=>psps
 this%mpi_enreg=>mpi_enreg
 this%crystal=>crystal

 this%wvl=>wvl
 this%results_gs=>results_gs

 this%pawtab=>pawtab
 this%pawrad=>pawrad
 this%pawang=>pawang
 this%pawfgr=>pawfgr
 this%npwarr=>npwarr
 this%cprj=>cprj

 this%occ=>occ
 this%rhor=>rhor
 this%rhog=>rhog
 this%pawrhoij=>pawrhoij

 this%nfftf=>nfftf
 this%mcg=>mcg
 this%mcprj=>mcprj
 this%my_natom=>my_natom
 this%cg=>cg
 this%kg=>kg
 this%ylm=>ylm
 this%ylmgr=>ylmgr

 this%ecore=>ecore

 this%dens_mat_real=>dens_mat_real
 if(present(emb_pot)) then
   this%has_embpot = .true.
   this%emb_pot=>emb_pot
 else
   this%has_embpot = .false.
 endif

 if(present(fock_mat)) then
   this%save_fock_mat = .true.
   this%fock_mat => fock_mat
 else
   this%save_fock_mat = .false.
 endif

 if(present(subham_sub))then
   this%subham_sub => subham_sub
 else
   ABI_ALLOCATE(this%subham_sub,(dim_all,dim_all))
   this%subham_sub = czero
   do ii=dim_sub+1,dim_all
     this%subham_sub(ii,ii) = cone
   enddo
 endif

 if(present(ireadwf)) then
   this%ireadwf = ireadwf
 else
   this%ireadwf = 0 !default
 endif


!useless
 this%taug=>taug
 this%taur=>taur
 this%paw_dmft=>paw_dmft
 this%dtefield=>dtefield
 this%electronpositron=>electronpositron

 this%pwind=>pwind
 this%pwind_alloc=>pwind_alloc
 this%pwnsfac=>pwnsfac

end subroutine subscf_init



!!****f* m_subscf/subscf_destroy
!! NAME
!! subscf_destroy
!!
!! FUNCTION
!! destroy subscf object
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
subroutine subscf_destroy(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_destroy'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout):: this

 this%dtfil=>null()
 this%dtset=>null()
 this%psps=>null()
 this%mpi_enreg=>null()
 this%crystal=>null()

 this%wvl=>null()
 this%results_gs=>null()

 this%pawtab=>null()
 this%pawrad=>null()
 this%pawang=>null()
 this%pawfgr=>null()
 this%npwarr=>null()

 this%occ=>null()

 this%nfftf=>null()
 this%mcg=>null()
 this%cg=>null()
 this%kg=>null()
 this%ylm=>null()
 this%ylmgr=>null()

 this%ecore=>null()

 if(this%has_embpot) this%emb_pot=>null()

 this%electronpositron=>null()

 this%taur=>null()
 this%taug=>null()

end subroutine subscf_destroy



!!****f* m_subscf/subscf_run
!! NAME
!! subscf_run
!!
!! FUNCTION
!! subscf driver
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
subroutine subscf_run(this,can2sub,dim_can,dim_sub,dim_all,hdr,&
&                     prtden,crystal_tot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_run'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout):: this
 type(hdr_type),intent(inout) :: hdr

 integer, intent(in) :: dim_can, dim_sub, dim_all
 integer, intent(in), optional :: prtden
 type(crystal_t),intent(in),optional ::crystal_tot

 complex(dpc), intent(in) :: can2sub(dim_can,dim_all)

 integer :: initialized,prtden_local

 initialized = 0

 if(present(prtden))then
   prtden_local = prtden
 else
   prtden_local = -1
 endif

 if(present(crystal_tot))then
   call subscf_core(this,this%dtset,this%crystal,this%psps,this%pawtab,this%pawrad,this%pawang,this%pawfgr,&
&   this%mpi_enreg,initialized,this%cprj,this%mcprj,&
&   this%nfftf,this%ecore,this%rhog,this%rhor,can2sub,dim_can,dim_sub,dim_all,hdr,prtden_local,crystal_tot)
 else
   call subscf_core(this,this%dtset,this%crystal,this%psps,this%pawtab,this%pawrad,this%pawang,this%pawfgr,&
&   this%mpi_enreg,initialized,this%cprj,this%mcprj,&
&   this%nfftf,this%ecore,this%rhog,this%rhor,can2sub,dim_can,dim_sub,dim_all,hdr,prtden_local)
 endif

end subroutine subscf_run


!!****f* m_subscf/subscf_core
!! NAME
!! subscf_core
!!
!! FUNCTION
!! main routine for scf calculation in subspace
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
subroutine subscf_core(this,dtset,crystal,psps,pawtab,pawrad,pawang,pawfgr,mpi_enreg,initialized,&
&                      cprj,mcprj,nfftf,ecore,rhog,rhor,can2sub,dim_can,dim_sub,dim_all,&
&                      hdr,prtden,crystal_tot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_core'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout) :: this
 type(dataset_type),intent(inout) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 type(crystal_t),intent(inout) :: crystal
 type(pseudopotential_type),intent(in) :: psps

 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawang_type), intent(in) :: pawang
 type(pawfgr_type), intent(inout) :: pawfgr
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)

 type(hdr_type),intent(inout) :: hdr
 
 integer,intent(inout) :: initialized,nfftf,mcprj
 type(pawcprj_type),pointer,intent(inout) :: cprj(:,:)

 real(dp),intent(inout) :: rhog(2,nfftf)
 real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)

 integer, intent(in) :: dim_can,dim_sub,dim_all,prtden
 type(crystal_t),intent(in),optional ::crystal_tot

 complex(dpc), intent(in) :: can2sub(dim_can,dim_all)

 real(dp),intent(in) :: ecore

 integer,parameter :: response=0
 integer :: conv_retcode,istep_updatedfock,mcg_sub
 integer :: ii,jj,iscf10,ispden,itypat,sz1,sz2
 integer :: nfftot,nfftotf,quit,prtfor,prtxml
 integer :: npwdiel,afford,dielstrt,choice,computed_forces,errid
 integer :: initialized0,istep_mix,ispmix,nfftmix,nfftmix_per_nfft,npawmix
 integer :: cplex,cplex_hf,mgfftf,ipositron,optene,optres,ngrvdw,nkxc,dbl_nnsclo,denpot
 integer :: has_dijhat,has_dijnd,has_dijfock,has_dijU,has_vhartree
 integer :: iatom,istep,nstep,n1xccc,n3xccc,nzlmopt,option,optxc
 integer :: nhatgrdim,ider,idir,izero,ipert,moved_atm_inside,moved_rhor
 integer :: my_natom,usexcnhat,usefock,usecprj,forces_needed,stress_needed,use_hybcomp
 integer :: optcut,optgr0,optgr1,optgr2,optrad
 integer :: mband_cprj,my_nspinor,ctocprj_choice,iorder_cprj
 integer :: wfmixalg,optcut_hf,optgr0_hf,optgr1_hf,optgr2_hf,optrad_hf,spare_mem,history_size,usecg,istep_fock_outer
 type(scf_history_type) :: scf_history_wf
 real(dp) :: boxcut,compch_fft,compch_sph,gsqcut,ecut,ecutf,ucvol,ucvol_local
 real(dp) :: zion,vxcavg,hyb_mixing,hyb_mixing_sr
 real(dp) :: etotal,deltae,elast,diffor,maxfor
 real(dp) :: res2,residm

 logical :: dummy_nhatgr
 logical :: reset_mixing=.false.
 logical,save :: tfw_activated=.false.
 character(len=1500) :: message

 type(energies_type) :: energies 
 type(ab7_mixing_object) :: mix

 !arrays
 integer :: ngfft(18),ngfftf(18),ngfftmix(18)
 integer,allocatable :: kg_diel(:,:)
 integer,allocatable :: dimcprj(:),dimcprj_srt(:),indsym(:,:,:),symrec(:,:,:),irrzon(:,:,:)
 real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),ph1d(:,:),ph1df(:,:)
 integer,allocatable :: l_size_atm(:)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: favg(3),gmet(3,3),gprimd(3,3),rmet(3,3),strsxc(6),vpotzero(2),red_ptot(3),tollist(12)
 real(dp) :: vnew_mean(dtset%nspden), vres_mean(dtset%nspden)
 real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:),phnons(:,:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vxc(:,:),vxc_hybcomp(:,:),vxctau(:,:,:),workr(:,:),xccc3d(:),ylmdiel(:,:)
 real(dp),allocatable :: grchempottn(:,:),grewtn(:,:),grhf(:,:),grnl(:),grvdw(:,:),grxc(:,:)
 real(dp),allocatable :: cg_sub(:,:),cg_new(:,:)

 real(dp) :: dielar(7)
 real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:),nvresid(:,:),susmat(:,:,:,:,:),tauresid(:,:),synlgr(:,:)
 real(dp),allocatable :: resid(:)

 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)

 type(fock_type),pointer :: fock

 type(pawcprj_type),allocatable, target :: cprj_local(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)

 !debug
 integer,save :: icalled = 0
! type(crystal_t),intent(in) ::crystal_tot

 !useless
 real(dp) :: qphon(3),rhopsg_dummy(0,0),rhopsr_dummy(0,0),rhor_dummy(0,0)
 type(efield_type) :: dtefield
 real(dp),allocatable :: doccde(:)
 type(ebands_t) :: ebands

 write(std_out,*) "Entering subscf_core"
 icalled = icalled + 1


 conv_retcode = 0
 initialized0 = initialized
 quit = 0

 dbl_nnsclo = 0
 deltae=zero ; elast=zero 

 ABI_ALLOCATE(this%eig_sub,(dim_all))
 this%eig_sub=zero

 ipert=0;idir=0;cplex=1
 istep_mix=1
 ipositron = 0
 forces_needed = 0; stress_needed = 0
 moved_atm_inside = 0
 vpotzero(:)=zero

 nstep=dtset%nstep
 my_natom = mpi_enreg%my_natom
 ecut=dtset%ecut
 ecutf=ecut; if (psps%usepaw==1) ecutf=dtset%pawecutdg
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg

 zion=zero
 do iatom=1,dtset%natom
   zion=zion+psps%ziontypat(dtset%typat(iatom))
 end do


 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
   mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
   mgfftf=dtset%mgfft;ngfftf(:)=ngfft(:)
 end if

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 nfftotf=product(ngfftf(1:3))

 call metric(gmet,gprimd,-1,rmet,crystal%rprimd,ucvol)
 ucvol_local = ucvol

 usecprj=0
 if (mcprj>0) then
  usecprj=1
 end if

 iscf10=mod(dtset%iscf,10)
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff
 tollist(8)=dtset%vdw_df_threshold

 optres=merge(0,1,dtset%iscf<10)

 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
 ABI_ALLOCATE(vhartr,(nfftf))
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vpsp,(nfftf))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vxctau,(nfftf,dtset%nspden,4*dtset%usekden))

 wfmixalg=dtset%fockoptmix/100
 use_hybcomp=0
 if(mod(dtset%fockoptmix,100)==11)use_hybcomp=1
 ABI_ALLOCATE(vxc_hybcomp,(nfftf,dtset%nspden*use_hybcomp))

 call energies_init(energies)
 select case(dtset%usepotzero)
 case(0,1)
   energies%e_corepsp   = ecore / ucvol
   energies%e_corepspdc = zero
 case(2)
   energies%e_corepsp   = zero
   energies%e_corepspdc = zero
 end select


 usefock=dtset%usefock
 nullify(fock)

 strsxc=zero
 ABI_ALLOCATE(grchempottn,(3,dtset%natom))
 grchempottn(:,:)=zero
 ABI_ALLOCATE(grewtn,(3,dtset%natom))
 ngrvdw=0
 ABI_ALLOCATE(grvdw,(3,ngrvdw))

 ABI_ALLOCATE(grnl,(3*dtset%natom))
 ABI_ALLOCATE(grxc,(3,dtset%natom))
 ABI_ALLOCATE(synlgr,(3,dtset%natom))

 nkxc=0
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=2*min(dtset%nspden,2)-1
 if (nkxc>0) then
   MSG_ERROR('nkxc>0, NYI')
   call check_kxc(dtset%ixc,dtset%optdriver)
 end if
 ABI_ALLOCATE(kxc,(nfftf,nkxc))

 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))


!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (dtset%iscf>=0) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
   dielar(7)=dtset%diemix;if (dtset%iscf>=10) dielar(7)=dtset%diemixmag
   ABI_ALLOCATE(nvresid,(nfftf,dtset%nspden))
   ABI_ALLOCATE(tauresid,(nfftf,dtset%nspden*dtset%usekden))
   if (nstep==0) then
     nvresid=zero
     tauresid=zero
   end if
   ABI_ALLOCATE(dtn_pc,(3,dtset%natom))
!  The next arrays are needed if iscf==5 and ionmov==4,
!  but for the time being, they are always allocated
   ABI_ALLOCATE(grhf,(3,dtset%natom))
!  Additional allocation for mixing within PAW
   npawmix=0
   if(psps%usepaw==1) then
     do iatom=1,my_natom
       itypat=this%pawrhoij(iatom)%itypat
       this%pawrhoij(iatom)%use_rhoijres=1
       sz1=this%pawrhoij(iatom)%cplex*pawtab(itypat)%lmn2_size
       sz2=this%pawrhoij(iatom)%nspden
       ABI_ALLOCATE(this%pawrhoij(iatom)%rhoijres,(sz1,sz2))
       do ispden=1,this%pawrhoij(iatom)%nspden
         this%pawrhoij(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(this%pawrhoij(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       this%pawrhoij(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       this%pawrhoij(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+this%pawrhoij(iatom)%nspden*pawtab(itypat)%lmnmix_sz*this%pawrhoij(iatom)%cplex
     end do
   end if
   if (dtset%iscf > 0) then
     denpot = AB7_MIXING_POTENTIAL
     if (dtset%iscf > 10) denpot = AB7_MIXING_DENSITY
     if (psps%usepaw==1.and.dtset%pawmixdg==0 .and. dtset%usewvl==0) then
       ispmix=AB7_MIXING_FOURRIER_SPACE;nfftmix=dtset%nfft;ngfftmix(:)=ngfft(:)
     else
       ispmix=AB7_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
     end if
     !TRangel: added to avoid segfaults with Wavelets
     nfftmix_per_nfft=0;if(nfftf>0) nfftmix_per_nfft=(1-nfftmix/nfftf)
     call ab7_mixing_new(mix, iscf10, denpot, ispmix, nfftmix, dtset%nspden, npawmix, errid, message, dtset%npulayit)
     if (errid /= AB7_NO_ERROR) then
       MSG_ERROR(message)
     end if
     if (dtset%mffmem == 0) then
       call ab7_mixing_use_disk_cache(mix, this%dtfil%fnametmp_fft)
     end if
!   else if (dtset%iscf==0.and.dtset%usewvl==1) then
!     ispmix=AB7_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
   end if
 else
   ABI_ALLOCATE(nvresid,(0,0))
   ABI_ALLOCATE(tauresid,(0,0))
   ABI_ALLOCATE(dtn_pc,(0,0))
   ABI_ALLOCATE(grhf,(0,0))
 end if ! iscf>0

 dielstrt = 0

 npwdiel = 1 
 afford = 0
 ABI_ALLOCATE(dielinv,(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
 ABI_ALLOCATE(susmat,(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
 ABI_ALLOCATE(kg_diel,(3,npwdiel))

 computed_forces=0
 choice=1 ; diffor=zero ; res2=zero
 etotal = zero
 !useless
 maxfor=zero
 ABI_ALLOCATE(fcart,(3,dtset%natom))
 ABI_ALLOCATE(forold,(3,dtset%natom))
 ABI_ALLOCATE(fred,(3,dtset%natom))
 ABI_ALLOCATE(gresid,(3,dtset%natom))

 fcart(:,:) = zero; forold(:,:) = zero; fred(:,:) = zero; gresid(:,:)=zero
 !end useless

 ABI_ALLOCATE(resid,(dtset%mband*dtset%nkpt*dtset%nsppol)) 
 resid = zero
 residm = zero

 prtfor=0;prtxml=0
 call scprqt(choice,dtset%cpus,deltae,diffor,dtset,&
& this%eig_sub,etotal,favg,fcart,energies%e_fermie,this%dtfil%fnameabo_app_eig,&
& this%dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
& maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
& this%occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
& psps%usepaw,vxcavg,dtset%wtk,crystal%xred,conv_retcode)


 !symmetry related stuffs
 ABI_ALLOCATE(irrzon,(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon(:,:,:)=0
 phnons(:,:,:)=zero
 indsym(:,:,:)=0
 symrec(:,:,:)=0
 if (dtset%nsym>1) then
   call setsym(indsym,irrzon,dtset%iscf,dtset%natom,&
&   nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,crystal%xred)

!  Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)
 else
!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1
 end if

!is this the right place?
 call symmetrize_xred(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,crystal%xred)
!Get cut-off for g-vectors
 if (psps%usepaw==1) then
   call wrtout(std_out,' FFT (fine) grid used in SCF cycle:','COLL')
 end if
!Compute large sphere cut-off (gsqcut):
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)
!Compute structure factor phases:
 call getph(crystal%atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,crystal%xred)
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(crystal%atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,crystal%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if


 if(psps%usepaw==1) then
   ABI_ALLOCATE(nhat,(nfftf,dtset%nspden*psps%usepaw))
   if (nstep==0) nhat=zero
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(my_natom))
   if (my_natom>0) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,dtset%typat,&
&      mpi_atmtab=mpi_enreg%my_atmtab)
     call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat,&
&      mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     ABI_DEALLOCATE(l_size_atm)
   end if
   compch_fft=-1.d5
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   write(std_out,*) "usexcnhat = ",usexcnhat
   if (usexcnhat==0.and.dtset%ionmov==4.and.dtset%iscf<10) then
     message = 'You cannot simultaneously use ionmov=4 and such a PAW psp file!'
     MSG_ERROR(message)
   end if

   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   call paw_an_nullify(paw_an)
   call paw_ij_nullify(paw_ij)
   has_dijhat=0;if (dtset%iscf==22) has_dijhat=1
   has_vhartree=0; if (dtset%prtvha > 0 .or. dtset%prtvclmb > 0) has_vhartree=1
   has_dijfock=0; if (usefock==1) has_dijfock=1
   has_dijnd=0;if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
   has_dijU=0; if (dtset%usepawu==5.or.dtset%usepawu==6) has_dijU=1
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,0,0,dtset%nspden,&
&   cplex,dtset%pawxcdev,dtset%typat,pawang,pawtab,has_vxc=1,has_vxc_ex=1,has_vhartree=has_vhartree,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,&
&   dtset%pawspnorb,dtset%natom,dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijfock=has_dijfock,has_dijhartree=1,has_dijnd=has_dijnd,has_dijso=1,has_dijhat=has_dijhat,&
&   has_dijU=has_dijU,has_pawu_occ=1,has_exexch_pot=1,nucdipmom=dtset%nucdipmom,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   compch_sph=-1.d5

!   ABI_ALLOCATE(dimcprj,(dtset%natom))
   ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
!   call pawcprj_getdim(dimcprj    ,dtset%natom,crystal%nattyp,dtset%ntypat,dtset%typat,pawtab,'R')
   call pawcprj_getdim(dimcprj_srt,dtset%natom,crystal%nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
   do itypat=1,dtset%ntypat
     if (pawtab(itypat)%usepawu>0) MSG_ERROR('PAW+U NYI!')!lpawumax=max(pawtab(itypat)%lpawu,lpawumax)
   end do
!   if (dtset%usedmatpu/=0.and.lpawumax>0) then
!     if (2*lpawumax+1/=size(dmatpawu,1).or.2*lpawumax+1/=size(dmatpawu,2)) then
!       message = 'Incorrect size for dmatpawu!'
!       MSG_BUG(message)
!     end if
!   end if

!  Allocation of projected WF (optional)
   if (usecprj==1) then
     iorder_cprj=0
     if (usefock==1) then
       ctocprj_choice = 1
       if (dtset%optforces == 1) then
        ctocprj_choice = 2; ! ncpgr = 3 
       end if
     end if
   endif

   !nullify(pawrhoij_ep);if(associated(electronpositron))pawrhoij_ep=>electronpositron%pawrhoij_ep
   !nullify(lmselect_ep);if(associated(electronpositron))lmselect_ep=>electronpositron%lmselect_ep

 else
!   ABI_ALLOCATE(dimcprj,(0))
   ABI_ALLOCATE(dimcprj_srt,(0))
   ABI_ALLOCATE(nhat,(0,0))
   ABI_DATATYPE_ALLOCATE(paw_ij,(0))
   ABI_DATATYPE_ALLOCATE(paw_an,(0))
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(0))
 endif !end usepaw


!copy cg to cg_sub
 mcg_sub = this%npwarr(1)*dim_sub !gamma point only
 ABI_ALLOCATE(cg_sub,(2,mcg_sub))
 call cg_zcopy(mcg_sub,this%cg,cg_sub)


 if(this%ireadwf==1)then
   !new cg
   ABI_ALLOCATE(cg_new,(2,mcg_sub))
   call cgtosub(cg_new,cg_sub,this%npwarr(1),this%subham_sub(1:dim_sub,1:dim_sub),dim_sub,dim_sub) !FIXME 
   call cg_zcopy(mcg_sub,cg_new,this%cg)
   ABI_DEALLOCATE(cg_new)

   !new rho
   if (psps%usepaw==1) then
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     call mkrho(this%cg,dtset,gprimd,irrzon,this%kg,this%mcg,&
&      mpi_enreg,this%npwarr,this%occ,this%paw_dmft,phnons,rhowfg,rhowfr,crystal%rprimd,1,ucvol,this%wvl%den,this%wvl%wfs)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
   else
     call mkrho(this%cg,dtset,gprimd,irrzon,this%kg,this%mcg,&
&      mpi_enreg,this%npwarr,this%occ,this%paw_dmft,phnons,rhog,rhor,crystal%rprimd,1,ucvol,this%wvl%den,this%wvl%wfs)
   end if

   if (psps%usepaw==1) then
     my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
     mband_cprj=dtset%mband
     if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
     mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
     ABI_DATATYPE_ALLOCATE(cprj_local,(dtset%natom,mcprj))
     ctocprj_choice = 1
     call pawcprj_alloc(cprj_local,0,dimcprj_srt)
     cprj=> cprj_local
     iatom=0 ; iorder_cprj=0 !0: sorted; 1:un-sorted
     idir = 0
     call ctocprj(crystal%atindx,this%cg,ctocprj_choice,cprj_local,gmet,gprimd,&
&     iatom,idir,iorder_cprj,dtset%istwfk,this%kg,dtset%kptns,&
&     this%mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&     dtset%mpw,dtset%natom,crystal%nattyp,dtset%nband,dtset%natom,ngfft,&
&     dtset%nkpt,dtset%nloalg,this%npwarr,dtset%nspinor,dtset%nsppol,&
&     dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&     ucvol,this%dtfil%unpaw,crystal%xred,this%ylm,this%ylmgr)

     pawrhoij_unsym=>this%pawrhoij
 
     call pawmkrhoij(crystal%atindx,crystal%atindx1,cprj,dimcprj_srt,dtset%istwfk,dtset%kptopt,&
&     dtset%mband,mband_cprj,mcprj,dtset%mkmem,mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,&
&     dtset%nspinor,dtset%nsppol,this%occ,dtset%paral_kgb,this%paw_dmft,dtset%pawprtvol,pawrhoij_unsym,&
&     this%dtfil%unpaw,dtset%usewvl,dtset%wtk)

     ABI_DEALLOCATE(dimcprj_srt)
!     ABI_DEALLOCATE(dimcprj)
     call pawcprj_free(cprj_local)
     ABI_DATATYPE_DEALLOCATE(cprj_local)
   end if
 endif


 istep_updatedfock=0

 !scf iteration
 do istep=1,max(1,nstep)

   if (istep==1) then
       !call symmetrize_xred(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,crystal%xred)

!      Get cut-off for g-vectors
!       if (psps%usepaw==1) then
!         call wrtout(std_out,' FFT (fine) grid used in SCF cycle:','COLL')
!       end if
!       call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!      Compute structure factor phases and large sphere cut-off (gsqcut):
!       call getph(crystal%atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,crystal%xred)

!       if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
!         call getph(crystal%atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,crystal%xred)
!       else
!         ph1df(:,:)=ph1d(:,:)
!       end if

     if (psps%usepaw==1) then
!      Check for non-overlapping spheres
       call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,crystal%xred)

!      Identify parts of the rectangular grid where the density has to be calculated
       optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
       if (forces_needed==1.or.(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)) then
         optgr1=dtset%pawstgylm
         if (stress_needed==1)  optrad=1
         if (dtset%pawprtwf==1) optrad=1
       end if

       call nhatgrid(crystal%atindx1,gmet,my_natom,dtset%natom,&
&        crystal%nattyp,ngfftf,psps%ntypat,optcut,optgr0,optgr1,optgr2,optrad,&
&        pawfgrtab,pawtab,crystal%rprimd,dtset%typat,ucvol,crystal%xred,&
&        comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&        comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)

       if(this%ireadwf==1)then
!        Compute n_tild + n_hat
         qphon(:)=zero
         call pawmkrho(1,compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&          my_natom,dtset%natom,dtset%nspden,dtset%nsym,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&          dtset%pawprtvol,this%pawrhoij,pawrhoij_unsym,pawtab,qphon,rhowfg,rhowfr,rhor,crystal%rprimd,dtset%symafm,&
&          symrec,dtset%typat,ucvol,dtset%usewvl,crystal%xred,rhog=rhog,pawnhat=nhat)
         ABI_DEALLOCATE(rhowfg)
         ABI_DEALLOCATE(rhowfr)
       endif
     end if
     !write(std_out,*) "debug: compch_fft = ", compch_fft
     !write(std_out,*) "debug: maxval(nhat) = ", maxval(nhat)
     
     moved_rhor = 0

   end if !end istep==1

   !Initialize/Update data in the case of an Exact-exchange (Hartree-Fock) or hybrid XC calculation
   hyb_mixing=zero;hyb_mixing_sr=zero
   if (usefock==1) then
     if (istep==1) then
       ! Initialize data_type fock for the calculation
       cplex_hf=cplex
       if (psps%usepaw==1) cplex_hf=dtset%pawcpxocc
       call fock_init(crystal%atindx,cplex_hf,dtset,fock,gsqcut,this%kg,mpi_enreg,crystal%nattyp,this%npwarr,pawang,pawfgr,pawtab,crystal%rprimd)
       if (fock%fock_common%usepaw==1) then
         optcut_hf = 0 ! use rpaw to construct local_pawfgrtab
         optgr0_hf = 0; optgr1_hf = 0; optgr2_hf = 0 ! dont need gY terms locally
         optrad_hf = 1 ! do store r-R
         call nhatgrid(crystal%atindx1,gmet,dtset%natom,dtset%natom,crystal%nattyp,ngfftf,psps%ntypat,&
&         optcut_hf,optgr0_hf,optgr1_hf,optgr2_hf,optrad_hf,fock%fock_common%pawfgrtab,pawtab,&
&         crystal%rprimd,dtset%typat,ucvol,crystal%xred,typord=1)
         iatom=-1;idir=0
         call ctocprj(crystal%atindx,this%cg,ctocprj_choice,cprj,gmet,gprimd,iatom,idir,&
&         iorder_cprj,dtset%istwfk,this%kg,dtset%kptns,this%mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,crystal%nattyp,dtset%nband,dtset%natom,ngfft, dtset%nkpt,dtset%nloalg,this%npwarr,dtset%nspinor,&
&         dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,this%dtfil%unpaw,&
&         crystal%xred,this%ylm,this%ylmgr)
       end if
       if(wfmixalg/=0)then
         spare_mem=0
         if(spare_mem==1)history_size=wfmixalg ! Not yet coded
         if(spare_mem==0)history_size=2*(wfmixalg-1)+1
!        Specific case of simple mixing : always history_size=1
         if(wfmixalg==2)history_size=1
         scf_history_wf%history_size=history_size
         usecg=2
         call scf_history_init(dtset,mpi_enreg,usecg,scf_history_wf)
       end if
     end if

     !Fock energy
     energies%e_exactX=zero
     if (fock%fock_common%optfor) then
       fock%fock_common%forces=zero
     end if

     if (istep==1 .or. istep_updatedfock==fock%fock_common%nnsclo_hf) then

       istep_updatedfock=1

       !Possibly mix the wavefunctions from different steps before computing the Fock operator
       if(wfmixalg/=0 .and. .not. (wfmixalg==2 .and. abs(scf_history_wf%alpha-one)<tol8) )then
         call wf_mixing(crystal%atindx1,this%cg,cprj,dtset,istep_fock_outer,this%mcg,mcprj,mpi_enreg,&
&         crystal%nattyp,this%npwarr,pawtab,scf_history_wf,0)
         istep_fock_outer=istep_fock_outer+1
       endif

       ! Update data relative to the occupied states in fock
       call fock_updatecwaveocc(this%cg,cprj,dtset,fock,indsym,this%mcg,mcprj,mpi_enreg,crystal%nattyp,this%npwarr,this%occ,ucvol)
       ! Possibly (re)compute the ACE operator
       if(fock%fock_common%use_ACE/=0) then
         call fock2ACE(this%cg,cprj,fock,dtset%istwfk,this%kg,dtset%kptns,dtset%mband,this%mcg,mcprj,dtset%mgfft,&
&         dtset%mkmem,mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,dtset%natom,dtset%nband,dtset%nfft,ngfft,dtset%nkpt,dtset%nloalg,this%npwarr,dtset%nspden,&
&         dtset%nspinor,dtset%nsppol,dtset%ntypat,this%occ,dtset%optforces,paw_ij,pawtab,ph1d,psps,crystal%rprimd,&
&         dtset%typat,usecprj,dtset%use_gpu_cuda,dtset%wtk,crystal%xred,this%ylm)
       end if

       !Should place a test on whether there should be the final exit of the istep loop.
       !This test should use focktoldfe.
       !This should update the info in fock%fock_common%fock_converged.
       !For the time being, fock%fock_common%fock_converged=.false. , so the loop end with the maximal value of nstep always,
       !except when nnsclo_hf==1 (so the Fock operator is always updated), in which case, the usual exit tests (toldfe, tolvrs, etc)
       !work fine.
       !if(fock%fock_common%nnsclo_hf==1 .and. fock%fock_common%use_ACE==0)then
       if(fock%fock_common%nnsclo_hf==1)then
         fock%fock_common%fock_converged=.TRUE.
       end if


       !Depending on fockoptmix, possibly restart the mixing procedure for the potential
       if(mod(dtset%fockoptmix,10)==1)then
         istep_mix=1
       end if
     else
       istep_updatedfock=istep_updatedfock+1
     end if

     !Used locally
     hyb_mixing=fock%fock_common%hyb_mixing ; hyb_mixing_sr=fock%fock_common%hyb_mixing_sr

   end if ! usefock


   if (istep==1.or.(mod(dtset%fockoptmix,100)==11 .and. istep_updatedfock==1)) then
!    PAW only: we sometimes have to compute compensation density
!    and eventually add it to density from WFs
     nhatgrdim=0
     dummy_nhatgr = .False.
     if (psps%usepaw==1.and.this%ireadwf==0.and.usexcnhat==0) then
       nhatgrdim=0;if (dtset%xclevel==2) nhatgrdim=usexcnhat*dtset%pawnhatxc
       ider=2*nhatgrdim;izero=0
       if (nhatgrdim>0)   then
         ABI_ALLOCATE(nhatgr,(cplex*nfftf,dtset%nspden,3*nhatgrdim))
       else
         ABI_ALLOCATE(nhatgr,(0,0,0))
         dummy_nhatgr = .True.
       end if
       call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
&       nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,&
&       nhatgr,nhat,this%pawrhoij,this%pawrhoij,pawtab,k0,crystal%rprimd,ucvol_local,dtset%usewvl,crystal%xred,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&       distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
!       if (this%dtfil%ireadwf/=0.and.this%dtfil%ireadden==0) then
!         rhor(:,:)=rhor(:,:)+nhat(:,:)
!         call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
!       end if
     end if

     optene = 4 * optres
     if(dtset%iscf==-3) optene=4

     if (.not.allocated(nhatgr))  then
       ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3*nhatgrdim))
       dummy_nhatgr = .True.
     end if

     call setvtr(crystal%atindx1,dtset,energies,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcut,&
&     istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&     crystal%nattyp,nfftf,ngfftf,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,&
&     n1xccc,n3xccc,optene,pawrad,pawtab,ph1df,psps,rhog,rhor,rmet,crystal%rprimd,&
&     strsxc,ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,this%wvl,&
&     xccc3d,crystal%xred,electronpositron=this%electronpositron,&
&     taug=this%taug,taur=this%taur,vxc_hybcomp=vxc_hybcomp,vxctau=vxctau,add_tfw=tfw_activated)

     ! set the zero of the potentials here
     if(dtset%usepotzero==2) then
       vpsp(:) = vpsp(:) + ecore / ( zion * ucvol )
     end if

     if ((nhatgrdim>0.and.nstep>0).or.dummy_nhatgr) then
       ABI_DEALLOCATE(nhatgr)
     end if

!debug
!     if(icalled == 3) then 
!     if(present(hdr)) then
!       call fftdatar_write("vpsp",this%dtfil%fnameabo_app_vpsp,dtset%iomode,hdr,&
!&        crystal_tot,ngfftf,1,nfftf,dtset%nspden,vpsp,mpi_enreg)
!     endif
!    else
!       bstruct_tmp = ebands_from_dtset(dtset, this%npwarr)
!       call hdr_init(bstruct_tmp,codvsn,dtset,hdr_tmp,pawtab,0,psps,this%wvl%descr,&
!&        comm_atom=this%mpi_enreg%comm_atom,mpi_atmtab=this%mpi_enreg%my_atmtab)
!       call ebands_free(bstruct_tmp)
!       call fftdatar_write("vpsp",this%dtfil%fnameabo_app_vpsp,dtset%iomode,hdr_tmp,&
!&        crystal,ngfftf,1,nfftf,dtset%nspden,vpsp,mpi_enreg)
!     endif
!     endif
!end debug
   endif


!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------
!  PAW: Compute energies and potentials in the augmentation regions (spheres)
!  Compute pseudopotential strengths (Dij quantities)
   if (psps%usepaw==1)then

!    Local exact exch.: impose occ. matrix if required
     if (dtset%useexexch>0) then
       MSG_ERROR('useexexch not supported!')
!       call setrhoijpbe0(dtset,initialized0,istep,istep_mix,&
!&       spaceComm,my_natom,dtset%natom,dtset%ntypat,pawrhoij,pawtab,dtset%typat,&
!&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if

!    Computation of on-site densities/potentials/energies
     nzlmopt=0;if (istep_mix==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep_mix>2) nzlmopt=dtset%pawnzlm
     call paw_an_reset_flags(paw_an) ! Force the recomputation of on-site potentials
     call paw_ij_reset_flags(paw_ij,self_consistent=.true.) ! Force the recomputation of Dij
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,dtset%ixc,my_natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,&
&     this%pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&     electronpositron=this%electronpositron,vpotzero=vpotzero)

     !write(std_out,*) "debug: compch_sph = ",compch_sph
     !do iatom=1,my_natom
     !   write(std_out,*) pawrhoij(iatom)%rhoijp(:,:)
     !enddo

!    Correct the average potential with the calculated constant vpotzero
!    Correct the total energies accordingly
!    vpotzero(1) = -beta/ucvol
!    vpotzero(2) = -1/ucvol sum_ij rho_ij gamma_ij
     write(message,'(a,f14.6,2x,f14.6)') &
&     ' average electrostatic smooth potential [Ha] , [eV]',SUM(vpotzero(:)),SUM(vpotzero(:))*Ha_eV
     call wrtout(std_out,message,'COLL')
     vtrial(:,:)=vtrial(:,:)+SUM(vpotzero(:))
     if(option/=1)then
!      Fix the direct total energy (non-zero only for charged systems)
       energies%e_paw=energies%e_paw-SUM(vpotzero(:))*dtset%charge
!      Fix the double counting total energy accordingly (for both charged AND
!      neutral systems)
       energies%e_pawdc=energies%e_pawdc-SUM(vpotzero(:))*zion+vpotzero(2)*dtset%charge
     end if


     if (dtset%usepawu>0.and.(ipositron/=1)) then
       MSG_ERROR('NYI')
     endif

     call pawdij(cplex,dtset%enunit,gprimd,ipert,my_natom,dtset%natom,nfftf,nfftotf,&
&     dtset%nspden,psps%ntypat,paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,&
&     pawrad,this%pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,k0,dtset%spnorbscl,&
&     ucvol_local,dtset%charge,vtrial,vxc,crystal%xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     mpi_comm_grid=mpi_enreg%comm_fft,&
&     hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&     nucdipmom=dtset%nucdipmom)

!    Symetrize Dij
     call symdij(gprimd,indsym,ipert,my_natom,dtset%natom,dtset%nsym,&
&     psps%ntypat,0,paw_ij,pawang,dtset%pawprtvol,pawtab,crystal%rprimd,dtset%symafm,symrec,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     if (has_dijhat==1) then
       call symdij(gprimd,indsym,ipert,my_natom,dtset%natom,dtset%nsym,&
&       psps%ntypat,1,paw_ij,pawang,dtset%pawprtvol,pawtab,crystal%rprimd,dtset%symafm,symrec,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if

   endif

   if (dtset%usepawu>0.and.dtset%macro_uj>0.and.istep>1.and.ipositron/=1) then
     MSG_ERROR('NYI!')
   endif

   if(nstep==0)exit

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   if(dtset%iscf>=0)then
     write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,message,'COLL')
   end if

   call subscf_vtorho(this,dtset,psps,crystal,mpi_enreg,this%dtfil,istep,compch_fft,&
&    cg_sub,mcg_sub,cprj,mcprj,usecprj,pawtab,pawfgr,pawfgrtab,pawang,paw_ij,this%pawrhoij,&
&    rhog,rhor,nhat,nvresid,optres,res2,tauresid,&
&    this%kg,this%ylm,this%ylmgr,vtrial,energies,ph1d,fock,my_natom,dtset%natom,psps%ntypat,0,nfftf,&
&    gmet,gprimd,indsym,symrec,irrzon,phnons,rmet,ucvol,this%paw_dmft,this%wvl,can2sub,dim_can,dim_sub,dim_all)

   if (dtset%iscf<0) exit

   if (dtset%iscf>=10) then
     optene = 1  ! use double counting scheme (default)
     if (dtset%iscf==22) optene = -1

!    Add the Fock contribution to E_xc and E_xcdc if required
     if (usefock==1) then
       energies%e_fockdc=two*energies%e_fock
     end if

!    if the mixing is the ODA mixing, compute energy and new density here
     if (dtset%iscf==22) then
        MSG_ERROR('iscf==22 is not supported')
!       call odamix(deltae,dtset,&
!&       elast,energies,etotal,gprimd,gsqcut,kxc,mpi_enreg,&
!&       my_natom,nfftf,ngfftf,nhat,nkxc,psps%ntypat,nvresid,n3xccc,optres,&
!&       paw_ij,paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
!&       red_ptot,psps,rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,&
!&       usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred,&
!&       taug=taug,taur=taur,vxctau=vxctau,add_tfw=tfw_activated)
     end if
!    If the density mixing is required, compute the total energy here
! TODO: add taur taug tauresid if needed

     !no orbital energy available for frozen orbitals
     if(dim_sub<dim_all) optene = 0
     call etotfor(crystal%atindx1,deltae,diffor,dtefield,dtset,&
&     elast,this%electronpositron,energies,&
&     etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,gresid,grewtn,grhf,grnl,grvdw,&
&     grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,my_natom,&
&     crystal%nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,psps%ntypat,nvresid,n1xccc,n3xccc,&
&     optene,computed_forces,optres,pawang,pawfgrtab,pawrad,this%pawrhoij,pawtab,&
&     ph1df,red_ptot,psps,rhog,rhor,rmet,crystal%rprimd,symrec,synlgr,ucvol,&
&     psps%usepaw,vhartr,vpsp,vxc,this%wvl%descr,this%wvl%den,xccc3d,crystal%xred)
   endif

   if (dtset%iscf>=10) then
!    Check exit criteria
     choice=2
     call scprqt(choice,dtset%cpus,deltae,diffor,dtset,&
&     this%eig_sub,etotal,favg,fcart,energies%e_fermie,this%dtfil%fnameabo_app_eig,&
&     this%dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     this%occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,crystal%xred,conv_retcode,&
&     electronpositron=this%electronpositron,fock=fock)

!    Check if we need to exit the loop
     if (istep==nstep) quit=1
!     quit_sum=quit
!     call xmpi_sum(quit_sum,spaceComm,ierr)
!     if (quit_sum>0) quit=1

!    If criteria in scprqt say to quit, then exit the loop over istep.
     if (quit==1) exit
   end if

   if (dtset%iscf>=10 .and.dtset%iscf/=22) then
     call newrho(crystal%atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,&
&     dtset,etotal,fcart,pawfgr%fintocoa,&
&     gmet,grhf,gsqcut,initialized,ispmix,istep_mix,kg_diel,kxc,&
&     mgfftf,mix,pawfgr%coatofin,moved_atm_inside,mpi_enreg,my_natom,crystal%nattyp,nfftf,&
&     nfftmix,nfftmix_per_nfft,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,nvresid,psps%ntypat,&
&     n1xccc,this%pawrhoij,pawtab,ph1df,psps,rhog,rhor,&
&     crystal%rprimd,susmat,psps%usepaw,vtrial,this%wvl%descr,this%wvl%den,crystal%xred,&
&     taug=this%taug,taur=this%taur,tauresid=tauresid)
   end if 

   optxc = 1

   if (dtset%iscf/=22) then
!    PAW: eventually recompute compensation density (and gradients)
     nhatgrdim=0
     if ( allocated(nhatgr) ) then
       ABI_DEALLOCATE(nhatgr)
     end if
     if (psps%usepaw==1) then
       ider=-1;if (dtset%iscf>=10.and.((dtset%xclevel==2.and.dtset%pawnhatxc>0).or.usexcnhat==0)) ider=0
       if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) ider=ider+2
       if (ipositron==1) ider=-1
       if (ider>0) then
         nhatgrdim=1
         ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3))
       else
         ABI_ALLOCATE(nhatgr,(0,0,0))
       end if
       if (ider>=0) then
         izero=0
         call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,nfftf,ngfftf,&
&         nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,nhatgr,nhat,&
&         this%pawrhoij,this%pawrhoij,pawtab,k0,crystal%rprimd,ucvol_local,dtset%usewvl,crystal%xred,&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&         comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&         distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
       end if
     else
       ABI_ALLOCATE(nhatgr,(0,0,0))
     end if
!    Compute new potential from the trial density

     optene=2*optres;if(psps%usepaw==1) optene=2

! TODO: check if tauresid is needed here too for potential residual in the future for MGGA potential mixing
     call rhotov(dtset,energies,gprimd,gsqcut,istep,kxc,mpi_enreg,nfftf,ngfftf, &
&     nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,optene,optres,optxc,&
&     rhog,rhor,crystal%rprimd,strsxc,ucvol_local,psps%usepaw,usexcnhat,&
&     vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,this%wvl,xccc3d,crystal%xred,&
&     electronpositron=this%electronpositron,taug=this%taug,taur=this%taur,vxctau=vxctau,&
&     vxc_hybcomp=vxc_hybcomp,add_tfw=tfw_activated)

   end if


!   if(VERBOSE)then
!     call wrtout(std_out,'*. END MINIMIZATION ITERATIONS',"COLL")
!   end if

   initialized = 1


   istep_mix=istep_mix+1
   if (reset_mixing) then
     istep_mix=1;reset_mixing=.false.
   end if

 enddo !istep

 if (dtset%iscf > 0) then
   call ab7_mixing_deallocate(mix)
 end if

 if (usefock==1)then
   if(wfmixalg/=0)then
     call scf_history_free(scf_history_wf)
   end if
 end if

 call cleanup(this%results_gs,energies,etotal)
 call hdr_update(hdr,hdr%bantot,etotal,energies%e_fermie,&
&  residm,crystal%rprimd,this%occ,this%pawrhoij,crystal%xred,dtset%amu_orig(:,1),&
&  comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!write density to file
 if(prtden>0) then
   ABI_MALLOC(doccde,(dim_all))
   doccde=zero
   call ebands_init(dim_all,ebands,dtset%nelect,doccde,this%eig_sub,hdr%istwfk,hdr%kptns,hdr%nband,&
&   hdr%nkpt,hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk,&
&   hdr%charge, hdr%kptopt, hdr%kptrlatt_orig, hdr%nshiftk_orig, hdr%shiftk_orig, &
&   hdr%kptrlatt, hdr%nshiftk, hdr%shiftk)

   ABI_FREE(doccde)

   ebands%fermie  = energies%e_fermie
   ebands%entropy = energies%entropy

   select case(prtden)
     case(1)
       call fftdatar_write("density",this%dtfil%fnameabo_impden,dtset%iomode,hdr,&
         crystal_tot,ngfftf,1,nfftf,dtset%nspden,rhor,mpi_enreg,ebands=ebands)
     case(2)
       call fftdatar_write("density",this%dtfil%fnameabo_bathden,dtset%iomode,hdr,&
         crystal_tot,ngfftf,1,nfftf,dtset%nspden,rhor,mpi_enreg,ebands=ebands)
   end select

   call ebands_free(ebands)
 endif


 if (psps%usepaw==1) then
   if (dtset%iscf>0) then
     do iatom=1,my_natom
       this%pawrhoij(iatom)%lmnmix_sz=0
       this%pawrhoij(iatom)%use_rhoijres=0
       ABI_DEALLOCATE(this%pawrhoij(iatom)%kpawmix)
       ABI_DEALLOCATE(this%pawrhoij(iatom)%rhoijres)
     end do
   end if
!   if (recompute_cprj) then
!     usecprj=0;mcprj=0
!     call pawcprj_free(cprj)
!     ABI_DATATYPE_DEALLOCATE(cprj_local)
!   end if
   call paw_an_free(paw_an)
   call paw_ij_free(paw_ij)
   call pawfgrtab_free(pawfgrtab)

   ABI_DATATYPE_DEALLOCATE(pawfgrtab)
   ABI_DATATYPE_DEALLOCATE(paw_an)
   ABI_DATATYPE_DEALLOCATE(paw_ij)
   ABI_DEALLOCATE(nhat)
 end if



end subroutine subscf_core



subroutine cleanup(results_gs,energies,etotal)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cleanup'
!End of the abilint section

 type(results_gs_type),intent(inout) :: results_gs
 type(energies_type),intent(in) :: energies
 real(dp),intent(in) :: etotal

 call energies_copy(energies,results_gs%energies)
 results_gs%etotal     =etotal


end subroutine cleanup


!!****f* m_subscf/subscf_vtorho
!! NAME
!! subscf_vtorho
!!
!! FUNCTION
!! potential to density
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
subroutine subscf_vtorho(this,dtset,psps,crystal,mpi_enreg,dtfil,istep,compch_fft,&
& cg_sub,mcg_sub,cprj,mcprj,usecprj,pawtab,pawfgr,pawfgrtab,pawang,paw_ij,pawrhoij,&
& rhog,rhor,nhat,nvresid,optres,nres2,tauresid,&
& kg,ylm,ylmgr,vtrial,energies,ph1d,fock,my_natom,natom,ntypat,optforces,nfftf,&
& gmet,gprimd,indsym,symrec,irrzon,phnons,rmet,ucvol,paw_dmft,wvl,can2sub,dim_can,dim_sub,dim_all)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_vtorho'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout) :: this
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(crystal_t),intent(in) :: crystal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(energies_type), intent(inout) :: energies
 type(fock_type),pointer, intent(inout) :: fock
 type(datafiles_type), intent(in) :: dtfil

 integer, intent(in) :: istep,my_natom,natom,ntypat,nfftf,optforces,optres,mcg_sub
 integer, intent(in) :: usecprj,mcprj
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden), compch_fft, nres2
 real(dp), intent(in) :: cg_sub(2,mcg_sub)
 real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp), intent(out) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp), intent(out) :: nvresid(nfftf,dtset%nspden),tauresid(nfftf,dtset%nspden*dtset%usekden)
 type(pawcprj_type),pointer,intent(inout) :: cprj(:,:)

 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
 type(pawang_type), intent(in) :: pawang
 type(pawfgr_type), intent(in) :: pawfgr
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)

 integer, intent(in) :: dim_can,dim_sub,dim_all
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: indsym(4,dtset%nsym,natom)
 real(dp), intent(in) :: ucvol
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 complex(dpc), intent(in) :: can2sub(dim_can,dim_all)
 type(paw_dmft_type),intent(inout) :: paw_dmft
 type(wvl_data), intent(inout) :: wvl

 integer :: mcg_new,ierr,mbdkpsp,iorder_cprj
 integer :: usecprj_local,istwf_k,cplex,ipert
 integer :: isppol,ikg,ilm,nkpg,dimffnl,ider,idir
 integer :: ikpt_loc,ikpt,nband_k,my_ikpt
 integer :: n1,n2,n3,n4,n5,n6,npw_k
 integer :: iband,ii,iscf
 integer :: mband_cprj,mcprj_local,my_nspinor,usefock_ACE !,mcprj_tmp
 integer,parameter :: tim_mkrho=2
 logical :: paral_atom,usefock
 real(dp) :: ar
 real(dp) :: nelect_frozen

 !arrays
 integer :: nband_sub(1)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: rhodum(1),kpoint(3),ylmgr_dum(0,0,0),qpt(3)
 real(dp), allocatable :: ylm_k(:,:),kinpw(:),kpg_k(:,:),subham(:)
 type(gs_hamiltonian_type) :: gs_hamk
 real(dp),allocatable :: cgrvtrial(:,:),vlocal(:,:,:,:),ffnl(:,:,:,:),ph3d(:,:,:),zshift(:)
 real(dp),allocatable :: cg_new(:,:), tmp_real(:,:),tmp_img(:,:)
 real(dp),allocatable :: rhowfg(:,:),rhowfr(:,:),focknk(:),ghc_dummy(:,:)
 complex(dpc),allocatable :: dens_mat(:,:),tmp(:,:)

! type(pawcprj_type),allocatable :: cprj_tmp(:,:)
 type(pawcprj_type),allocatable,target :: cprj_local(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:) 

 integer :: me_distrb,mpi_comm_sphgrid,nproc_distrb,spaceComm_distrb
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()

 !useless
 real(dp),allocatable :: doccde(:)

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 usefock = (dtset%usefock==1 .and. associated(fock))

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 paral_atom=(my_natom/=natom)
 if(paral_atom) MSG_ERROR('paral_atom not supported')
 compch_fft=-1.d5

 usefock = (dtset%usefock==1 .and. associated(fock))
 usefock_ACE=0
 if (usefock) usefock_ACE=fock%fock_common%use_ACE


 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   MSG_BUG('wrong values for nfft, nfftf!')
 end if

 if(optforces.ne.0) then
   MSG_BUG('optforces!=0')
 endif

 spaceComm_distrb=mpi_enreg%comm_cell
 if (mpi_enreg%paral_kgb==1) spaceComm_distrb=mpi_enreg%comm_kpt
 if (mpi_enreg%paral_hf ==1) spaceComm_distrb=mpi_enreg%comm_kpt
 nproc_distrb=xmpi_comm_size(spaceComm_distrb)
 me_distrb=xmpi_comm_rank(spaceComm_distrb)
 mpi_comm_sphgrid=mpi_enreg%comm_fft

 usecprj_local=0;if (psps%usepaw==1) usecprj_local=1

 iscf = dtset%iscf
! if(.not. wvlbigdft) then
   energies%e_eigenvalues = zero
   energies%e_kinetic     = zero
   energies%e_nonlocalpsp = zero
   if (usefock) then
     energies%e_fock=zero
     energies%e_fockdc=zero
   end if
!   grnl(:)=zero
!   resid(:) = zero ! JWZ 13 May 2010. resid and eigen need to be fully zeroed each time before use
!   eigen(:) = zero
!   bdtot_index=0
!   ibg=0;icg=0
   mbdkpsp=dtset%mband*dtset%nkpt*dtset%nsppol
! end if

   if (usefock) then
     ABI_ALLOCATE(focknk,(mbdkpsp))
     focknk=zero
   end if


 if(iscf>=0 .or. iscf==-3) then
   if (optres==1) then
     nvresid=rhor
     tauresid=this%taur
   end if
!  NC and plane waves
   if (psps%usepaw==0 .and. dtset%usewvl==0) then
     rhor=zero
!    PAW
   elseif(psps%usepaw==1) then
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     rhowfr(:,:)=zero
   end if
 end if

 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,natom,&
& dtset%typat,crystal%xred,dtset%nfft,dtset%mgfft,dtset%ngfft,crystal%rprimd,dtset%nloalg,&
& paw_ij=paw_ij,ph1d=ph1d,usecprj=usecprj_local,electronpositron=this%electronpositron,fock=fock,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 mcprj_local=0 ; mband_cprj=0
 if (psps%usepaw==1) then
   mband_cprj=dtset%mband
   if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
   iorder_cprj=0 ; mcprj_local=mcprj
   if (usecprj==0) then
     mcprj_local=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
     !This is a check but should always be true since scfcv allocated cprj
     !anyway
     if (allocated(cprj_local)) then
       !Was allocated in scfcv so we just destroy and reconstruct it as desired
       call pawcprj_free(cprj_local)
       ABI_DATATYPE_DEALLOCATE(cprj_local)
     end if
     ABI_DATATYPE_ALLOCATE(cprj_local,(dtset%natom,mcprj_local))
     call pawcprj_alloc(cprj_local,0,gs_hamk%dimcprj)
     cprj=> cprj_local
   end if
 end if


 ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamk%nvloc))


 do isppol=1,dtset%nsppol
   ikpt_loc = 0
   ikg = 0

   ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
   call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
   call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
   ABI_DEALLOCATE(cgrvtrial)

   call load_spin_hamiltonian(gs_hamk,isppol,vlocal=vlocal,with_nonlocal=.true.)

   ikpt = 0
   do while (ikpt_loc < dtset%nkpt)

     ikpt_loc = ikpt_loc + 1
     ikpt = ikpt_loc
     my_ikpt = mpi_enreg%my_kpttab(ikpt)

     !nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)

     npw_k=this%npwarr(ikpt)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,dim_sub,isppol,me_distrb)) then
!       eigen(1+bdtot_index : nband_k+bdtot_index) = zero
!       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     if (mpi_enreg%paral_kgb==1) my_bandfft_kpt => bandfft_kpt(my_ikpt)
     call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

     istwf_k=dtset%istwfk(ikpt)
     kpoint(:)=dtset%kptns(:,ikpt)

!     ABI_ALLOCATE(zshift,(nband_k))
!     zshift(:)=dtset%eshift 

     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Set up remaining of the Hamiltonian

!    Compute (1/2) (2 Pi)**2 (k+G)**2:
     ABI_ALLOCATE(kinpw,(npw_k))
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=3*optforces*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;idir=0;dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
       call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&         npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&         psps%usepaw,psps%useylm,ylm_k,ylmgr)
     end if


!    Load k-dependent part in the Hamiltonian datastructure
!       - Compute 3D phase factors
!       - Prepare various tabs in case of band-FFT parallelism
!       - Load k-dependent quantities in the Hamiltonian
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))

     if(usefock_ACE/=0) then
       call load_k_hamiltonian(gs_hamk,kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,&
&         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,fockACE_k=fock%fockACE(ikpt,isppol),ph3d_k=ph3d,&
&         compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1),&
&         compute_gbound=(mpi_enreg%paral_kgb/=1))
     else
       call load_k_hamiltonian(gs_hamk,kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,&
&         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,&
&         compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1),&
&         compute_gbound=(mpi_enreg%paral_kgb/=1))
     endif

!    Load band-FFT tabs (transposed k-dependent arrays)
     if (mpi_enreg%paral_kgb==1) then
         if (istep<=1) then
           call prep_bandfft_tabs(gs_hamk,ikpt,dtset%mkmem,mpi_enreg)
         end if
         call load_k_hamiltonian(gs_hamk,npw_fft_k=my_bandfft_kpt%ndatarecv, &
&         gbound_k =my_bandfft_kpt%gbound, &
&         kinpw_k  =my_bandfft_kpt%kinpw_gather, &
&         kg_k     =my_bandfft_kpt%kg_k_gather, &
&         kpg_k    =my_bandfft_kpt%kpg_k_gather, &
          ffnl_k   =my_bandfft_kpt%ffnl_gather, &
          ph3d_k   =my_bandfft_kpt%ph3d_gather)
     end if


     ! Setup gemm_nonlop
     if (gemm_nonlop_use_gemm) then
         !set the global variable indicating to gemm_nonlop where to get its data from
         gemm_nonlop_ikpt_this_proc_being_treated = my_ikpt
         if (istep <= 1) then
           !Init the arrays
           call make_gemm_nonlop(my_ikpt,gs_hamk%npw_fft_k,gs_hamk%lmnmax, &
&           gs_hamk%ntypat, gs_hamk%indlmn, gs_hamk%nattyp, gs_hamk%istwf_k, gs_hamk%ucvol, gs_hamk%ffnl_k,&
&           gs_hamk%ph3d_k)
         end if
     end if

     if (usefock) then
       call fock_updateikpt(fock%fock_common,ikpt,isppol)
     end if


     !mcg_sub = npw_k*dim_sub
     !ABI_ALLOCATE(cg_sub,(2,mcg_sub))
     !call cgtosub(cg_sub,this%cg,npw_k,can2sub,dim_can,dim_sub) !FIXME

     call subscf_mkham(this,dtset,mpi_enreg,gs_hamk,ikpt,my_nspinor,&
&      this%subham_sub(1:dim_sub,1:dim_sub),this%eig_sub(1:dim_sub),dim_sub,cg_sub,mcg_sub)


     !ABI_DEALLOCATE(cg_sub)

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ph3d)
!     ABI_DEALLOCATE(zshift)
     ABI_DEALLOCATE(kinpw)
   enddo
 enddo

 ABI_DEALLOCATE(vlocal)

 ABI_ALLOCATE(doccde,(dim_sub*dtset%nkpt*dtset%nsppol))
 doccde=zero

 call xmpi_bcast(this%eig_sub,0,mpi_enreg%comm_kpt,ierr)
! write(std_out,*) 'debug eigen:'
! write(std_out,*) this%eig_sub
 call xmpi_bcast(this%subham_sub,0,mpi_enreg%comm_kpt,ierr)

! call newocc(doccde,this%eig_sub,energies%entropy,energies%e_fermie,dtset%spinmagntarget,&
!& dtset%mband,dtset%nband,dtset%nelect,dtset%nkpt,dtset%nspinor,&
!& dtset%nsppol,this%occ,dtset%occopt,dtset%prtvol,dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)

 nband_sub = dim_sub

 nelect_frozen = zero
 do ii=dim_sub+1,dim_all
   nelect_frozen = nelect_frozen + this%occ(ii)
 enddo

 if(dtset%occopt>=3.and.dtset%occopt<=8)then
 call newocc(doccde,this%eig_sub(1:dim_sub),energies%entropy,energies%e_fermie,dtset%spinmagntarget,&
& dim_sub,nband_sub,dtset%nelect-nelect_frozen,dtset%nkpt,dtset%nspinor,&
& dtset%nsppol,this%occ(1:dim_sub),dtset%occopt,dtset%prtvol,dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)
 endif

 ABI_DEALLOCATE(doccde)

 ABI_ALLOCATE(dens_mat,(dim_sub,dim_sub))
 dens_mat = czero
 do ii=1,dim_sub
   dens_mat(ii,ii) = cmplx(this%occ(ii),kind=dp)
 enddo

 ABI_ALLOCATE(tmp,(dim_sub,dim_sub))
 call zgemm('N','C',dim_sub,dim_sub,dim_sub,cone,dens_mat,dim_sub,this%subham_sub(1:dim_sub,1:dim_sub),dim_sub,czero,tmp,dim_sub)
 call zgemm('N','N',dim_sub,dim_sub,dim_sub,cone,this%subham_sub(1:dim_sub,1:dim_sub),dim_sub,tmp,dim_sub,czero,dens_mat,dim_sub)
 ABI_DEALLOCATE(tmp)

 ABI_ALLOCATE(tmp_img,(dim_sub,dim_sub))
 tmp_img = aimag(dens_mat)
! write(std_out,*) 'max dens_mat_img=',maxval(abs(tmp_img))
 if(maxval(abs(tmp_img))>tol8) then
   MSG_WARNING('Density matrix is not real!') 
   write(std_out,*) "max abs imaginary element:", maxval(abs(tmp_img))
!   do ii=1,dim_sub
!     write(std_out,'(*(F9.6))') tmp_img(ii,:)
!   enddo
 endif
 ABI_DEALLOCATE(tmp_img)

 this%dens_mat_real(1,:,:) = real(dens_mat,kind=dp)
 this%dens_mat_real(2,:,:) = zero
 if(maxval(abs(aimag(dens_mat))) > tol8) then
   this%dens_mat_real(2,:,:) = aimag(dens_mat)
 endif
 ABI_DEALLOCATE(dens_mat)

!get cg_new
! ABI_ALLOCATE(tmp,(dim_can,dim_all))
! call zgemm('N','N',dim_can,dim_all,dim_all,cone,can2sub,dim_can,this%subham_sub,dim_all,czero,tmp,dim_can)
 ikpt = 1
 npw_k=this%npwarr(ikpt)
! mcg_new = npw_k*dim_all
 ABI_ALLOCATE(cg_new,(2,mcg_sub))
 call cgtosub(cg_new,cg_sub,npw_k,this%subham_sub(1:dim_sub,1:dim_sub),dim_sub,dim_sub) !FIXME 
 call cg_zcopy(mcg_sub,cg_new,this%cg)
 ABI_DEALLOCATE(cg_new)
! if(dim_all > dim_sub) cg_new(:,1+npw_k*dim_sub:npw_k*dim_all) = this%cg(:,1+npw_k*dim_sub:npw_k*dim_all)
! ABI_DEALLOCATE(tmp) 
!end get cg_new

 ABI_ALLOCATE(kinpw,(npw_k))
 call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg,kinpw,kpoint,npw_k,0,0)
 !FIXME
 !compute kinetic energy
 if(iscf>0 .or. iscf==-3)then
   do iband=1,dim_all
     if(abs(this%occ(iband))>tol8) then
       call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&       this%cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),&
&       this%cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),0)

       energies%e_kinetic = energies%e_kinetic + this%occ(iband)*ar
       energies%e_eigenvalues = energies%e_eigenvalues + this%occ(iband)*this%eig_sub(iband)
     endif
   end do
 endif

!compute fock energy contribution
 if(usefock)then
   energies%e_fock = zero
   ABI_ALLOCATE(ghc_dummy,(2,npw_k*my_nspinor))
   call fock_updateikpt(fock%fock_common,1,1)
   do iband=1,dim_all
     call fock_set_ieigen(gs_hamk%fockcommon,iband)
     call fock_ACE_getghc(this%cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),ghc_dummy,gs_hamk,mpi_enreg)
   enddo
   ABI_DEALLOCATE(ghc_dummy)

   focknk(1:dim_all) = fock%fock_common%eigen_ikpt(1:dim_all)

   do iband=1,dim_all
     if(abs(this%occ(iband))>tol8)then
       energies%e_fock=energies%e_fock + half*focknk(iband)*this%occ(iband)
     end if
   enddo
   
   ABI_DEALLOCATE(focknk)
 endif



 call mkrho(this%cg,dtset,gprimd,irrzon,kg,this%mcg,mpi_enreg,this%npwarr,this%occ,paw_dmft,phnons,&
& rhowfg,rhowfr,crystal%rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)



 if (iscf>0.or.iscf==-3) then
!  PAW: Build new rhoij quantities from new occ then symetrize them
!  Compute and add the compensation density to rhowfr to get the total density
   if (psps%usepaw==1) then
     if (paral_atom) then
        MSG_ERROR('paral_atom not supported')
!       ABI_DATATYPE_ALLOCATE(pawrhoij_unsym,(natom))
!       nspden_rhoij=pawrhoij_get_nspden(dtset%nspden,dtset%nspinor,dtset%pawspnorb)
!       call pawrhoij_alloc(pawrhoij_unsym,dtset%pawcpxocc,nspden_rhoij,dtset%nspinor,&
!&       dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0)
     else
       pawrhoij_unsym => pawrhoij
     end if
     usecprj_local = 0 !have cg_new, need cprj
     if (usecprj_local==1) then
        MSG_BUG('should not be here!')
!       call pawmkrhoij(atindx,atindx1,cprj,gs_hamk%dimcprj,dtset%istwfk,dtset%kptopt,dtset%mband,mband_cprj,&
!&       mcprj_local,dtset%mkmem,mpi_enreg,natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
!&       occ,dtset%paral_kgb,paw_dmft,dtset%pawprtvol,pawrhoij_unsym,dtfil%unpaw,&
!&       dtset%usewvl,dtset%wtk)
     else
!       mcprj_tmp=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
!       ABI_DATATYPE_ALLOCATE(cprj_tmp,(natom,mcprj_tmp))
!       call pawcprj_alloc(cprj_tmp,0,gs_hamk%dimcprj)
       call ctocprj(crystal%atindx,this%cg,1,cprj,gmet,gprimd,0,0,0,dtset%istwfk,kg,dtset%kptns,&
&       this%mcg,mcprj_local,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&       dtset%natom,crystal%nattyp,dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg,&
&       this%npwarr,dtset%nspinor,dtset%nsppol,ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&       ucvol,dtfil%unpaw,crystal%xred,ylm,ylmgr_dum)
       call pawmkrhoij(crystal%atindx,crystal%atindx1,cprj,gs_hamk%dimcprj,dtset%istwfk,dtset%kptopt,&
&       dtset%mband,mband_cprj,mcprj_local,dtset%mkmem,mpi_enreg,natom,dtset%nband,dtset%nkpt,&
&       dtset%nspinor,dtset%nsppol,this%occ,dtset%paral_kgb,paw_dmft,dtset%pawprtvol,pawrhoij_unsym,&
&       dtfil%unpaw,dtset%usewvl,dtset%wtk)
!       call pawcprj_free(cprj_tmp)
!       ABI_DATATYPE_DEALLOCATE(cprj_tmp)
     end if
!    Build symetrized packed rhoij and compensated pseudo density
     cplex=1;ipert=0;idir=0;qpt(:)=zero
     if(dtset%usewvl==0) then
       call pawmkrho(1,compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&       my_natom,natom,dtset%nspden,dtset%nsym,ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&       dtset%pawprtvol,pawrhoij,pawrhoij_unsym,pawtab,qpt,rhowfg,rhowfr,rhor,crystal%rprimd,dtset%symafm,&
&       symrec,dtset%typat,ucvol,dtset%usewvl,crystal%xred,pawnhat=nhat,rhog=rhog)
     endif
   endif

!  Find and print minimum and maximum total electron density and locations
!  Compute density residual (if required) and its squared norm
   if (iscf>=0) then
     if (psps%usepaw==0) then
       call prtrhomxmn(std_out,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor,ucvol=ucvol)
     else
       call prtrhomxmn(std_out,mpi_enreg,nfftf,pawfgr%ngfft,dtset%nspden,1,rhor,ucvol=ucvol)
     end if
     if (optres==1) then
       nvresid=rhor-nvresid
       call sqnorm_v(1,nfftf,nres2,dtset%nspden,optres,nvresid,mpi_comm_sphgrid=mpi_enreg%comm_fft)
       tauresid=this%taur-tauresid
     end if
   end if

 endif

 if (psps%usepaw==1) then
   if (usecprj==0) then
     call pawcprj_free(cprj_local)
     ABI_DATATYPE_DEALLOCATE(cprj_local)
   end if
 end if


 if(psps%usepaw==1.and.(iscf>=0.or.iscf==-3))  then
   ABI_DEALLOCATE(rhowfr)
   ABI_DEALLOCATE(rhowfg)
 endif

 ABI_DEALLOCATE(kinpw)
! ABI_DEALLOCATE(cg_new)

 call destroy_hamiltonian(gs_hamk)

end subroutine subscf_vtorho


!!****f* m_subscf/subscf_mkham
!! NAME
!! subscf_mkham
!!
!! FUNCTION
!! build subspace hamiltonian
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
subroutine subscf_mkham(this,dtset,mpi_enreg,gs_hamk,ikpt,my_nspinor,subham_sub,eig_sub,nband_k,&
&                       cg,mcg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_mkham'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout):: this
 type(gs_hamiltonian_type), intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

 integer,intent(in) :: ikpt,my_nspinor,nband_k,mcg
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(inout) :: eig_sub(nband_k)
 complex(dpc),intent(inout) :: subham_sub(nband_k,nband_k)

!local variables
 integer,parameter :: icg = 0
 integer :: npw_k
 integer :: iband1,iband2,isubh
 integer :: lwork,info

 real(dp),allocatable :: subham(:)
 real(dp), allocatable :: rwork(:)
 complex(dpc), allocatable :: zwork(:)

 npw_k=this%npwarr(ikpt)

 ABI_ALLOCATE(subham,(nband_k*(nband_k+1)))
 subham(:) = zero
 call wftosubham(subham,gs_hamk,mpi_enreg,cg,mcg,icg,nband_k,dtset%nbdblock,npw_k,my_nspinor,dtset%prtvol)

 isubh = 1
 do iband1=1,nband_k
   do iband2=1,iband1
     if(iband2/=iband1)then
       !if(abs(subham(isubh+1))>tol8)then
       !  write(std_out,*) 'complex H:', subham(isubh+1)
       !endif
       subham_sub(iband2,iband1)=cmplx(subham(isubh),subham(isubh+1),kind=dp)
       subham_sub(iband1,iband2)=cmplx(subham(isubh),-subham(isubh+1),kind=dp)
     else
       if(abs(subham(isubh+1))>tol8) MSG_ERROR('Hamiltonian is not Hermitian!') 
       subham_sub(iband2,iband1)=cmplx(subham(isubh),kind=dp) 
     endif
     isubh=isubh+2
   enddo
 enddo

! do iband1=1,nband_k
!   do iband2=1,iband1
!     if(abs(subham_sub(iband2,iband1)-conjg(subham_sub(iband1,iband2)))>tol8 ) then
!       MSG_ERROR('Hamiltonian is not Hermitian!')
!     endif
!   enddo
! enddo

 if(this%save_fock_mat)then
   do iband1=1,nband_k
     do iband2=1,nband_k
       this%fock_mat(iband1,iband2) = subham_sub(iband1,iband2)
     enddo
   enddo

   if(maxval(abs(aimag(this%fock_mat))) < tol8) then
     do iband1=1,nband_k
       do iband2=1,nband_k
         this%fock_mat(iband1,iband2) = cmplx(real(this%fock_mat(iband1,iband2)),kind=dp)
       enddo
     enddo
   endif
 endif


 !write(std_out,*) 'fock_mat'
 !write(std_out,*) this%fock_mat

 if(this%has_embpot) then
   do iband1=1,nband_k
     do iband2=1,nband_k
       subham_sub(iband1,iband2) = subham_sub(iband1,iband2) &
&        +cmplx(this%emb_pot(1,iband1,iband2),this%emb_pot(2,iband1,iband2),kind=dp)
     enddo
   enddo
 endif

 ABI_ALLOCATE(rwork,(3*nband_k-2))
 lwork = 65*nband_k ! Value to optimize speed of the diagonalization
 ABI_ALLOCATE(zwork,(lwork))
 call zheev('v','u',nband_k,subham_sub,nband_k,eig_sub,zwork,lwork,rwork,info)
 if(info.ne.0) MSG_ERROR('Hamiltonian diagonalization failed!')

 write(std_out,*) "Eigenvalues in subspace DFT:"
 do iband1=1,nband_k
   write(std_out,'(f20.12)') eig_sub(iband1)
 enddo

 ABI_DEALLOCATE(subham)
 ABI_DEALLOCATE(rwork)
 ABI_DEALLOCATE(zwork)

end subroutine subscf_mkham


!!****f* m_subscf/wftosubham
!! NAME
!! wftosubham
!!
!! FUNCTION
!! build subspace hamiltonian
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
!! subscf_mkham
!!
!! SOURCE
subroutine wftosubham(subham,gs_hamk,mpi_enreg,cg,mcg,icg,nband,nbdblock,npw,nspinor,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wftosubham'
!End of the abilint section

 implicit none

 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(inout) :: mpi_enreg

 integer,intent(in) :: icg,mcg
 integer,intent(in) :: npw,nband,nbdblock,nspinor,prtvol

 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: subham(nband*(nband+1))
! complex(dpc),intent(out) :: subham(nband,nband)

!local variables
 integer,parameter :: cpopt=-1, tim_getghc=1, use_vnl=0, use_subovl=0, igsc=0
 integer :: iblock,nblock,ibandmin,ibandmax,iband
 integer :: icg_shift,sij_opt,istwf_k,isubh,isubo
 integer :: mgsc=1
 integer :: spacedim,blockdim,ioff,ii

 real(dp) :: eval,lambda_dum
 real(dp),allocatable :: cwavef(:,:), ghc(:,:), gvnlc(:,:)
 real(dp) :: gsc_dummy(2,1),subovl_dummy(0),subvnl_dummy(0)

 type(pawcprj_type) :: cprj_dum(1,1)

! type(xg_t) :: subsub
! type(xgBlock_t) :: X
! type(xgBlock_t) :: AX

 ABI_ALLOCATE(ghc,(2,npw*nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw*nspinor))

 istwf_k=gs_hamk%istwf_k

 isubh=1;isubo=1

 sij_opt = 0

! spacedim = npw*nspinor
! blockdim=mpi_enreg%nproc_band*mpi_enreg%bandpp 
! blockdim=1
! ABI_ALLOCATE(ghc,(2,spacedim*nband))
! ABI_ALLOCATE(ghc,(2,spacedim*blockdim))
! ABI_ALLOCATE(gvnlc,(2,spacedim*blockdim))

! write(std_out,*) 'blockdim=',blockdim
! write(std_out,*) 'nband=',nband

 print*,'nblock=',nblock
 print*,'nbdblock=',nbdblock
 print*,'nband=',nband
 print*,'me=',xmpi_comm_rank(xmpi_world)

 nblock=(nband-1)/nbdblock+1
! nblock = nband/blockdim
 ! Loop over blocks of bands. In the standard band-sequential algorithm, nblock=nband.
! ioff = 0
 do iblock=1,nblock

   ! Loop over bands in a block
   ! This loop can be MPI-parallelized, over processors attached to the same k point
   ibandmin=1+(iblock-1)*nbdblock
   ibandmax=min(iblock*nbdblock,nband)

!   if (mpi_enreg%paral_kgb==0) then
!     call multithreaded_getghc(cpopt,cg(:,1+ioff:blockdim*spacedim+ioff),cprj_dum,&
!&     ghc(:,1+ioff:blockdim*spacedim+ioff),gsc_dummy,gs_hamk,gvnlc,lambda_dum,mpi_enreg,blockdim,prtvol,sij_opt,tim_getghc,0)
!   else
!     call prep_getghc(cg(:,1+ioff:blockdim*spacedim+ioff),gs_hamk,gvnlc,&
!&     ghc(:,1+ioff:blockdim*spacedim+ioff),gsc_dummy,lambda_dum,blockdim,&
!&     mpi_enreg,prtvol,sij_opt,cpopt,cprj_dum,already_transposed=.false.)
!   endif

   ! Big iband loop
   do iband=ibandmin,ibandmax
     icg_shift=npw*nspinor*(iband-1)+icg

     call cg_zcopy(npw*nspinor,cg(1,1+icg_shift),cwavef)

!    By setting ieigen to iband, Fock contrib. of this iband to the energy will be calculated
     call fock_set_ieigen(gs_hamk%fockcommon,iband)

     call getghc(cpopt,cwavef,cprj_dum,ghc,gsc_dummy,gs_hamk,gvnlc,&
&     eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)
   enddo

   call mksubham(cg,ghc,gsc_dummy,gvnlc,iblock,icg,igsc,istwf_k,&
&   isubh,isubo,mcg,mgsc,nband,nbdblock,npw,&
&   nspinor,subham,subovl_dummy,subvnl_dummy,use_subovl,use_vnl,mpi_enreg%me_g0)

!   ioff = ioff + blockdim*spacedim
 enddo

 ABI_DEALLOCATE(gvnlc)

! call xg_init(subsub,SPACE_C,nband,nband) 
! call xgBlock_map(X,cg,SPACE_C,spacedim,nband)
! call xgBlock_map(AX,ghc,SPACE_C,spacedim,nband)
! call xgBlock_gemm(X%trans,AX%normal,1.0d0,X,AX,0.d0,subsub%self)

! call xg_get_cmplx_array(subsub,subham)
! call xg_free(subsub)

 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(cwavef)

end subroutine wftosubham


!!****f* m_subscf/cgtosub
!! NAME
!! cgtosub
!!
!! FUNCTION
!! transform cg to subspace
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
!!
!! SOURCE
subroutine cgtosub(cg_sub,cg,npw,can2sub,dim_can,dim_sub)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgtosub'
!End of the abilint section

 implicit none

 integer,intent(in) :: npw,dim_can,dim_sub
 real(dp),intent(in) :: cg(2,npw*dim_can)
 real(dp),intent(inout) :: cg_sub(2,npw*dim_sub)
 complex(dpc),intent(in) :: can2sub(dim_can,dim_sub)

!local variables
 real(dp),allocatable :: can2sub_re(:,:),can2sub_im(:,:) 


 ABI_ALLOCATE(can2sub_re,(dim_can,dim_sub))
 ABI_ALLOCATE(can2sub_im,(dim_can,dim_sub))
 can2sub_re = real(can2sub,kind=dp)
 can2sub_im = aimag(can2sub)


 call dgemm('N','N',npw,dim_sub,dim_can, one,cg(1,:),npw,can2sub_re,dim_can,zero,cg_sub(1,:),npw)
 call dgemm('N','N',npw,dim_sub,dim_can,-one,cg(2,:),npw,can2sub_im,dim_can, one,cg_sub(1,:),npw)
 call dgemm('N','N',npw,dim_sub,dim_can, one,cg(1,:),npw,can2sub_im,dim_can,zero,cg_sub(2,:),npw)
 call dgemm('N','N',npw,dim_sub,dim_can, one,cg(2,:),npw,can2sub_re,dim_can, one,cg_sub(2,:),npw)

 ABI_DEALLOCATE(can2sub_re)
 ABI_DEALLOCATE(can2sub_im)

end subroutine cgtosub


end module m_subscf
