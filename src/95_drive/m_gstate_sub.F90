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
 use libxc_functionals
 use m_errors
 use m_xmpi

 use m_gemm_nonlop
 use m_bandfft_kpt
 use m_efield

 use m_crystal,          only : crystal_t,crystal_init
 use defs_wvltypes,      only : wvl_data
 use m_geometry,         only : fixsym
 use m_kg,               only : kpgio, getph
 use m_berryphase_new,   only : init_e_field_vars
 use m_mkrho,            only : initro
 use m_initylmg,         only : initylmg
 use m_symtk,            only : matr3inv

 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_pawfgr,           only : pawfgr_type,pawfgr_init
 use m_pawrad,           only : pawrad_type
 use m_paw_occupancies,  only : initrhoij
 use m_pspini,           only : pspini
 use m_paw_init,         only : pawinit,paw_gencond
 use m_paw_sphharm,      only : setsym_ylm
 use m_spacepar,         only : setsym
 use m_common,           only : setup1

 use m_mpinfo,           only : proc_distrb_cycle

 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,paw_dmft_type
 use m_electronpositron, only : electronpositron_type,init_electronpositron,destroy_electronpositron

 use m_subscf

 implicit none

 private

 public :: gstate_sub

contains

subroutine gstate_sub(acell,dtset,psps,rprim,mpi_enreg,dtfil,wvl,&
& cg,pawtab,pawrad,pawang,xred,&
& dens,can2sub,dim_can,dim_sub,&
& emb_pot) !optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gstate_sub'
!End of the abilint section

 implicit none

 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(wvl_data),intent(inout) :: wvl !NYI



 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrad_type), intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type),intent(inout) :: pawang
 type(pawfgr_type) :: pawfgr


 integer, intent(in) :: dim_can, dim_sub

 real(dp),intent(inout) :: acell(3)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp), intent(in) :: cg(2,*)

 real(dp), intent(inout) :: dens(dim_sub,dim_sub),rprim(3,3)
 real(dp), intent(in), optional :: emb_pot(dim_sub,dim_sub)
 complex(dpc), intent(in) :: can2sub(dim_can,dim_sub)

 type(subscf_type) :: subscf_args
 type(crystal_t) :: crystal
 logical :: has_to_init,call_pawinit
 character(len=500) :: message

 type(paw_dmft_type) :: paw_dmft !NYI
 type(efield_type) :: dtefield !NYI
 type(electronpositron_type),pointer :: electronpositron !NYI
 !scalar
 integer :: comm,comm_psp,me,my_natom
 integer :: mgfftf,nfftf,nfftot
 integer :: bantot,psp_gencond,gnt_option,usecprj,option
 integer :: mcprj,mband_cprj,ncpgr,itypat
 integer :: cnt,spin,ikpt,band,my_nspinor,mcg
 integer :: use_sc_dmft,pwind_alloc !useless
 integer,parameter :: formeig=0,level=101,response=0 ,cplex1=1, master=0

 real(dp) :: ecut_eff,ecutdg_eff, gsqcut_eff,gsqcutc_eff,gsqcut_shp,hyb_range_fock
 real(dp) :: ecore,ucvol


 !arrays
 integer :: ngfft(18),ngfftf(18)
 integer,allocatable :: npwarr(:),kg(:,:),npwtot(:)
 integer,allocatable :: irrzon(:,:,:),indsym(:,:,:),symrec(:,:,:)
 integer,allocatable :: dimcprj_srt(:)
 integer,pointer :: pwind(:,:,:) !useless

 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_for_kg(3,3)
 real(dp) :: gmet(3,3),gprimd(3,3),gprimd_for_kg(3,3),gmet_for_kg(3,3)
 real(dp),allocatable :: occ(:),ph1df(:,:),ylm(:,:),ylmgr(:,:,:)
 real(dp),allocatable :: phnons(:,:,:)
 real(dp),pointer :: pwnsfac(:,:) !useless
 real(dp),pointer :: rhor(:,:),rhog(:,:),taug(:,:),taur(:,:)

 type(pawrhoij_type),pointer :: pawrhoij(:)

 DBG_ENTER("COLL")
!###########################################################
!### 01. Initializations XML, MPI, WVL, etc

!Init MPI data
 comm=mpi_enreg%comm_cell; me=xmpi_comm_rank(comm)

!Set up MPI information from the dataset
 my_natom=mpi_enreg%my_natom

!Define FFT grid(s) sizes (be careful !)
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

!###########################################################
!### 02. Calls setup1, kpgio, initylmg

 ecore=zero
 !results_gs%pel(1:3)   =zero
 !results_gs%grchempottn(:,:)=zero
 !results_gs%grewtn(:,:)=zero
 !call energies_init(results_gs%energies)

!Set up for iterations
 call setup1(acell,bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& dtset%natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,rprim,rprimd,ucvol,psps%usepaw)

!In some cases (e.g. getcell/=0), the plane wave vectors have
! to be generated from the original simulation cell
 rprimd_for_kg=rprimd
 !if (dtset%getcell/=0.and.dtset%usewvl==0) rprimd_for_kg=args_gs%rprimd_orig
 call matr3inv(rprimd_for_kg,gprimd_for_kg)
 gmet_for_kg=matmul(transpose(gprimd_for_kg),gprimd_for_kg)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 ABI_ALLOCATE(npwtot,(dtset%nkpt))
 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet_for_kg,dtset%istwfk,kg, &
&   dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
&   dtset%mpw,npwarr,npwtot,dtset%nsppol)
   call bandfft_kpt_init1(bandfft_kpt,dtset%istwfk,kg,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol)
 else
   ABI_ALLOCATE(kg,(0,0))
   npwarr = 0
   npwtot = 0
 end if

 !if(dtset%wfoptalg == 1 .and. psps%usepaw == 1) then
 !NYI
 !  call init_invovl(dtset%nkpt)  
 !end if

 if(dtset%use_gemm_nonlop == 1 .and. dtset%use_gpu_cuda/=1) then
   ! set global variable
   gemm_nonlop_use_gemm = .true.
   call init_gemm_nonlop(dtset%nkpt)
 else
   gemm_nonlop_use_gemm = .false.
 end if

!Set up the Ylm for each k point
 if ( dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   if (psps%useylm==1) then
     option=0
     if (dtset%prtstm==0.and.dtset%iscf>0.and.dtset%positron/=1) option=1 ! compute gradients of YLM
     if (dtset%berryopt==4 .and. dtset%optstress /= 0 .and. psps%usepaw==1) option = 1 ! compute gradients of YLM
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
&     psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
   end if
 else
   ABI_ALLOCATE(ylm,(0,0))
   ABI_ALLOCATE(ylmgr,(0,0,0))
 end if

 !FIXME
 has_to_init = .true.

!###########################################################
!### 03. Calls pspini

!Open and read pseudopotential files
 comm_psp=mpi_enreg%comm_cell
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,&
& pawrad,pawtab,psps,rprimd,comm_mpi=comm_psp)

 if (psps%usepaw==1.and.dtset%usefock==1)then
   do itypat=1,dtset%ntypat
     if(pawtab(itypat)%has_fock==0) then
       write(message, '(a)' )&
&       'The PAW data file does not contain Fock information. Change the PAW data file.'
       MSG_BUG(message)
     end if
   end do
 end if

!In case of isolated computations, ecore must set to zero
!because its contribution is counted in the ewald energy as the ion-ion interaction.
 if (dtset%icoulomb == 1) ecore = zero


!Initialize PAW atomic occupancies
 ABI_DATATYPE_ALLOCATE(pawrhoij,(my_natom*psps%usepaw))
 if (psps%usepaw==1.and.has_to_init) then
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,&
&   dtset%lpawu,my_natom,dtset%natom,dtset%nspden,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,pawrhoij,dtset%pawspnorb,pawtab,dtset%spinat,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if


 !simply use default occ
 ABI_ALLOCATE(occ,(dtset%mband*dtset%nkpt*dtset%nsppol))
 occ(:) = dtset%occ_orig(:,1)

!###########################################################
!### 04. Symmetry operations when nsym>1

!Do symmetry stuff only for nsym>1
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
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
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!  Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)
 else
!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1
 end if

!###########################################################
!### 05. Calls inwffil

 ! if paral_kgb == 0, it may happen that some processors are idle (no entry in proc_distrb)
 ! but mkmem == nkpt and this can cause integer overflow in mcg or allocation error.
 ! Here we count the number of states treated by the proc. if cnt == 0, mcg is then set to 0.
 cnt = 0
 do spin=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     do band=1,dtset%nband(ikpt + (spin-1) * dtset%nkpt)
       if (.not. proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, band, band, spin, mpi_enreg%me_kpt)) cnt = cnt + 1
     end do
   end do
 end do

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 mcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 if (cnt == 0) then
   mcg = 0
   write(message,"(2(a,i0))")"rank: ",mpi_enreg%me, "does not have wavefunctions to treat. Setting mcg to: ",mcg
   MSG_WARNING(message)
 end if

#if 0
 if (dtset%usewvl == 0 .and. dtset%mpw > 0 .and. cnt /= 0)then
   if (my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol > floor(real(HUGE(0))/real(dtset%mpw) )) then
     write (message,'(9a)')&
&     "Default integer is not wide enough to store the size of the wavefunction array (mcg).",ch10,&
&     "This usually happens when paral_kgb == 0 and there are not enough procs to distribute kpts and spins",ch10,&
&     "Action: if paral_kgb == 0, use nprocs = nkpt * nsppol to reduce the memory per node.",ch10,&
&     "If this does not solve the problem, use paral_kgb 1 with nprocs > nkpt * nsppol and use npfft/npband/npspinor",ch10,&
&     "to decrease the memory requirements. Consider also OpenMP threads."
     MSG_ERROR_NOSTOP(message,ii)
     write (message,'(5(a,i0), 2a)')&
&     "my_nspinor: ",my_nspinor, ", mpw: ",dtset%mpw, ", mband: ",dtset%mband,&
&     ", mkmem: ",dtset%mkmem, ", nsppol: ",dtset%nsppol,ch10,&
&     'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS...) compiled in int64 mode'
     MSG_ERROR(message)
   end if
 end if
#endif

 !ABI_ALLOCATE(eigen,(dtset%mband*dtset%nkpt*dtset%nsppol))
 !ABI_ALLOCATE(resid,(dtset%mband*dtset%nkpt*dtset%nsppol))
 !eigen(:)=zero ; resid(:)=zero


!Initialize wavefunctions.
!FIXME
!how to write wf to file?



!###########################################################
!### 07. Calls setup2

!Further setup
 !ABI_ALLOCATE(start,(3,dtset%natom))
 !call setup2(dtset,npwtot,start,wvl%wfs,xred)


!###########################################################
!### 08. Compute new occupation numbers
!FIXME

! results_gs%energies%entropy=zero

!###########################################################
!### 09. Generate an index table of atoms

!Definition of atindx array
!Generate an index table of atoms, in order for them to be used type after type.
! ABI_ALLOCATE(atindx,(dtset%natom))
! ABI_ALLOCATE(atindx1,(dtset%natom))
! ABI_ALLOCATE(nattyp,(psps%ntypat))
! indx=1
! do itypat=1,psps%ntypat
!   nattyp(itypat)=0
!   do iatom=1,dtset%natom
!     if(dtset%typat(iatom)==itypat)then
!       atindx(iatom)=indx
!       atindx1(indx)=iatom
!       indx=indx+1
!       nattyp(itypat)=nattyp(itypat)+1
!     end if
!   end do
! end do

 call crystal_init(dtset%amu_orig(:,1),crystal,dtset%spgroup,dtset%natom,dtset%npsp,&
& psps%ntypat,dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
& dtset%nspden==2.and.dtset%nsppol==1,.false.,psps%title,&
& symrel=dtset%symrel,tnons=dtset%tnons,symafm=dtset%symafm)


!Compute structure factor phases for current atomic pos:
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
 call getph(crystal%atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)

!###########################################################
!### 10. PAW related operations

!Initialize paw_dmft, even if neither dmft not paw are used
 use_sc_dmft=dtset%usedmft
 if(dtset%paral_kgb>0) use_sc_dmft=0
 call init_sc_dmft(dtset%nbandkss,dtset%dmftbandi,dtset%dmftbandf,dtset%dmft_read_occnd,dtset%mband,&
& dtset%nband,dtset%nkpt,dtset%nspden, &
& dtset%nspinor,dtset%nsppol,occ,dtset%usedmft,paw_dmft,use_sc_dmft,dtset%dmft_solv,mpi_enreg)

!PAW: 
!1- Initialize values for several arrays unchanged during iterations
!2- Initialize data for LDA+U
!3- Eventually open temporary storage file
 if(psps%usepaw==1) then
!  1-
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

!  Test if we have to call pawinit
!  Some gen-cond have to be added...
   call paw_gencond(dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1.or.call_pawinit) then
     gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     hyb_range_fock=zero;if (dtset%ixc<0) call libxc_functionals_get_hybridparams(hyb_range=hyb_range_fock)
     call pawinit(gnt_option,gsqcut_shp,hyb_range_fock,dtset%pawlcutd,dtset%pawlmix,&
&     psps%mpsang,dtset%pawnphi,dtset%nsym,dtset%pawntheta,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%xclevel,dtset%usepotzero)

     ! Update internal values
     call paw_gencond(dtset,gnt_option,"save",call_pawinit)
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   call setsym_ylm(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,rprimd,symrec,pawang%zarot)
!  2-Initialize and compute data for LDA+U, EXX, or LDA+DMFT
   pawtab(:)%usepawu=0
   pawtab(:)%useexexch=0
   pawtab(:)%exchmix=zero
 end if

!###########################################################
!### 11. Initialize (eventually) electron-positron data and
!###     electric and magnetic field data

!Initialize (eventually) electron-positron data
 nullify (electronpositron)
!###########################################################
! Initialisation of cprj
 usecprj=0; mcprj=0;mband_cprj=0
#if 0
 if (usecprj==0) then
   ABI_DATATYPE_ALLOCATE(cprj,(0,0))
 end if
 if (usecprj==1) then
   mband_cprj=dtset%mband;if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
   mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
!Was allocated above for valgrind sake so should always be true (safety)
   if (allocated(cprj)) then
     call pawcprj_free(cprj)
     ABI_DATATYPE_DEALLOCATE(cprj)
   end if
   ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
   ncpgr=0
   if (dtset%usefock==1) then
     if (dtset%optforces == 1) then
       ncpgr = 3
     end if
   end if
   ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
   call pawcprj_getdim(dimcprj_srt,dtset%natom,crystal%nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
   call pawcprj_alloc(cprj,ncpgr,dimcprj_srt)
 end if
#endif
!###########################################################
!### 12. Operations dependent of iscf value

!Get starting charge density : rhor as well as rhog
!Also initialize the kinetic energy density
 ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))
 ABI_ALLOCATE(taur,(nfftf,dtset%nspden*dtset%usekden))
 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(taug,(2,nfftf*dtset%usekden))

 !use atomic density
 if (has_to_init) then
   call initro(crystal%atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,&
&   mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,crystal%nattyp,nfftf,&
&   ngfftf,dtset%nspden,psps%ntypat,dtset%paral_kgb,psps,pawtab,ph1df,&
&   psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,psps%usepaw,&
&   dtset%ziontypat,dtset%znucl)
 end if
 ABI_DEALLOCATE(ph1df)
!###########################################################
!### 13. If needed, initialize SCF history variables

!If needed, initialize atomic density in SCF history

!Electric field: initialization stage
 call init_e_field_vars(dtefield,dtset,gmet,gprimd,kg,&
& mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&
& pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)



!  ###########################################################
!  ### 14. Run SCF

 if(present(emb_pot)) then
  call subscf_init(subscf_args,dtfil,dtset,psps,crystal,nfftf,&
&  pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&  ylm,ylmgr,kg,cg,mcg,my_natom,npwarr,ecore,wvl,&
&  occ,rhog,rhor,pawrhoij,dens,dim_sub,&
&  taug,taur,paw_dmft,dtefield,pwind_alloc,pwind,pwnsfac,electronpositron,&
&  emb_pot)
 else
  call subscf_init(subscf_args,dtfil,dtset,psps,crystal,nfftf,&
&  pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&  ylm,ylmgr,kg,cg,mcg,my_natom,npwarr,ecore,wvl,&
&  occ,rhog,rhor,pawrhoij,dens,dim_sub,&
&  taug,taur,paw_dmft,dtefield,pwind_alloc,pwind,pwnsfac,electronpositron)
 endif

 call subscf_run(subscf_args,can2sub,dim_can,dim_sub)

 call subscf_destroy(subscf_args)

end subroutine gstate_sub




end module m_gstate_sub
