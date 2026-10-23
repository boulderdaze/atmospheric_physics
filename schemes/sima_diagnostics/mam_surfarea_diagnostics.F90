! MAM aerosol surface area density (SAD) plumbing demonstration for the
! SAD -> MICM sandbox (NCAR/musica#1030).
!
! Calls the abstract aerosol interface (aero_state%surf_area_dens) for the
! MAM model the way CAM's aero_model_surfarea does for the tropospheric
! heterogeneous-chemistry rates in mo_usrrxt, and writes the results to
! history:
!   - SAD_TROP, REFF_TROP: total surface area density and effective radius
!     over the chemistry species types (CAM history names)
!   - sad_a<mm>, dmaer_a<mm>: per-mode surface area density and wet
!     number-mode diameter, the sfc/dm_aer pair CAM's hetrxtrate consumes
! It writes nothing into MICM: the SAD -> surface-reaction rate parameters are
! filled by mam_micm_surface_rates, which runs just before this scheme.
!
! It also checks that provider in place. Every MICM surface reaction carries
! two rate parameters per grid cell, '<label>.effective radius [m]' and
! '<label>.particle number concentration [# m-3]', and MICM's rate is
! 4 pi r^2 N / (r/D_g + 4/(v gamma)). CAM's hetrxtrate is the same expression
! per mode with 4 pi r^2 N replaced by the mode's surface area density and r
! by half its wet number-mode diameter. So for every surface reaction found in
! the mechanism this scheme recovers 4 pi r^2 N from the rate parameters,
! writes it to history (srx_sad_<label>, cm2/cm3), and writes the index of
! the MAM mode whose (sad_a<mm>, dmaer_a<mm>/2) pair it equals to within
! 1e-10 relative (srx_chk_<label>): 0 where both are zero (above the
! tropopause), -1 where nothing matches. Reactions mapped to other species
! types than CAM's chemistry list (sea salt) will not match by design.
module mam_surfarea_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: mam_surfarea_diagnostics_init
   public :: mam_surfarea_diagnostics_run

   ! Per-mode history field names, built at init.
   integer                        :: nmodes = 0
   character(len=16), allocatable :: sad_name(:)    ! sad_a<mm>
   character(len=16), allocatable :: dmaer_name(:)  ! dmaer_a<mm>

   ! MICM surface reactions found in the mechanism at init (rate-parameter
   ! labels '<label>.effective radius [m]' / '<label>.particle number
   ! concentration [# m-3]'), their rate-parameter slots and history names.
   real(kind_phys),  parameter :: match_rtol = 1.0e-10_kind_phys ! relative tolerance of the identity check
   integer                                          :: nsurf = 0
   character(len=:),   allocatable :: surf_label(:)
   integer,            allocatable :: surf_iradius(:)  ! slot of the effective radius [m]
   integer,            allocatable :: surf_inumber(:)  ! slot of the particle number concentration [# m-3]
   character(len=32),  allocatable :: srx_sad_name(:)  ! srx_sad_<label>: recovered 4 pi r^2 N (cm2/cm3)
   character(len=32),  allocatable :: srx_chk_name(:)  ! srx_chk_<label>: matching MAM mode index

contains

!> \section arg_table_mam_surfarea_diagnostics_init Argument Table
!! \htmlinclude mam_surfarea_diagnostics_init.html
   subroutine mam_surfarea_diagnostics_init(errmsg, errflg)
      use cam_history,       only: history_add_field
      use mam_mode_metadata, only: ntot_amode_val

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer          :: m
      character(len=3) :: trnum

      errmsg = ''
      errflg = 0

      nmodes = ntot_amode_val
      allocate(sad_name(nmodes), dmaer_name(nmodes), stat=errflg, errmsg=errmsg)
      if (errflg /= 0) then
         errmsg = 'mam_surfarea_diagnostics_init: allocate history names: '//trim(errmsg)
         return
      end if

      call history_add_field('SAD_TROP', &
           'tropospheric aerosol surface area density, chemistry species types', 'lev', 'avg', 'cm2/cm3')
      call history_add_field('REFF_TROP', &
           'tropospheric aerosol effective radius, chemistry species types', 'lev', 'avg', 'cm')

      ! Per-mode fields, CAM's per-mode naming pattern (write(trnum,'(i3.3)')m).
      do m = 1, nmodes
         write(trnum, '(i3.3)') m
         sad_name(m)   = 'sad_a'//trnum(2:3)
         dmaer_name(m) = 'dmaer_a'//trnum(2:3)

         call history_add_field(trim(sad_name(m)), &
              'aerosol surface area density, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'cm2/cm3')
         call history_add_field(trim(dmaer_name(m)), &
              'wet dgnum used for surface area density, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'cm')
      end do

      call find_micm_surface_reactions(errmsg, errflg)

   end subroutine mam_surfarea_diagnostics_init

   ! Enumerate the MICM surface reactions from the rate-parameter labels the
   ! mechanism configuration produced (available after musica_ccpp's register
   ! phase) and register their history fields. Field names are the label with
   ! any character outside [A-Za-z0-9_] replaced by '_', truncated to fit.
   subroutine find_micm_surface_reactions(errmsg, errflg)
      use cam_history,              only: history_add_field
      use cam_history_support,      only: fieldname_len
      use musica_ccpp_micm,         only: rate_parameters_ordering
      use musica_ccpp_micm_surface, only: find_surface_reactions_in_ordering

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: j, k, nlab
      character(len=fieldname_len)  :: stem

      call find_surface_reactions_in_ordering(rate_parameters_ordering, surf_label, &
           surf_iradius, surf_inumber, errmsg, errflg)
      if (errflg /= 0) return
      nsurf = size(surf_label)
      if (nsurf == 0) return

      allocate(srx_sad_name(nsurf), srx_chk_name(nsurf), stat=errflg, errmsg=errmsg)
      if (errflg /= 0) then
         errmsg = 'mam_surfarea_diagnostics_init: allocate surface reaction tables: '//trim(errmsg)
         return
      end if

      do j = 1, nsurf
         nlab = len_trim(surf_label(j))

         ! history names: 8-character prefix + sanitized label, within fieldname_len
         stem = ' '
         stem(1:min(nlab, fieldname_len - 8)) = sanitize(surf_label(j)(1:min(nlab, fieldname_len - 8)))
         srx_sad_name(j) = 'srx_sad_'//trim(stem)
         srx_chk_name(j) = 'srx_chk_'//trim(stem)
         do k = 1, j - 1
            if (srx_sad_name(k) == srx_sad_name(j)) then
               errflg = 1
               errmsg = 'mam_surfarea_diagnostics_init: MICM surface reaction labels "'// &
                        trim(surf_label(k))//'" and "'//trim(surf_label(j))// &
                        '" map to the same history field name '//trim(srx_sad_name(j))
               return
            end if
         end do

         call history_add_field(trim(srx_sad_name(j)), &
              'MICM surface reaction '//trim(surf_label(j))//': 4 pi r^2 N from its rate parameters', &
              'lev', 'avg', 'cm2/cm3')
         call history_add_field(trim(srx_chk_name(j)), &
              'MICM surface reaction '//trim(surf_label(j))// &
              ': MAM mode whose (SAD, wet dgnum/2) its (N, r) equal; 0 = both zero, -1 = mismatch', &
              'lev', 'avg', '1')
      end do

   end subroutine find_micm_surface_reactions

   pure function sanitize(str) result(out)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: out
      integer :: i
      out = str
      do i = 1, len(str)
         select case (str(i:i))
         case ('a':'z', 'A':'Z', '0':'9', '_')
         case default
            out(i:i) = '_'
         end select
      end do
   end function sanitize

!> \section arg_table_mam_surfarea_diagnostics_run Argument Table
!! \htmlinclude mam_surfarea_diagnostics_run.html
   subroutine mam_surfarea_diagnostics_run(ncol, pver, troplev, t, pmid, pi, rate_parameters, &
        sad_chem_spec_types, errmsg, errflg)
      use cam_history,            only: history_out_field
      use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                        aerosol_instances_get_state, &
                                        aerosol_instances_get_num_models
      use aerosol_properties_mod, only: aerosol_properties
      use aerosol_state_mod,      only: aerosol_state

      integer,          intent(in)  :: ncol
      integer,          intent(in)  :: pver
      integer,          intent(in)  :: troplev(:)  ! (ncol) tropopause vertical layer index
      real(kind_phys),  intent(in)  :: t(:,:)      ! (ncol,pver) air temperature [K]
      real(kind_phys),  intent(in)  :: pmid(:,:)   ! (ncol,pver) air pressure [Pa]
      real(kind_phys),  intent(in)  :: pi
      real(kind_phys),  intent(in)  :: rate_parameters(:,:,:)  ! (ncol,pver,nparams) MICM rate parameters
      character(len=*), intent(in)  :: sad_chem_spec_types(:)  ! aerosol species types making up the reactive surface (namelist)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: iaermod, m, j, i, k
      class(aerosol_properties), pointer :: aero_props
      class(aerosol_state),      pointer :: aero_state_obj

      integer         :: beglev(ncol), endlev(ncol)
      real(kind_phys) :: relhum(ncol, pver)
      real(kind_phys) :: sad(ncol, pver)             ! surface area density [cm2/cm3]
      real(kind_phys) :: reff(ncol, pver)            ! effective radius [cm]
      real(kind_phys) :: sfc(ncol, pver, nmodes)     ! per-mode surface area density [cm2/cm3]
      real(kind_phys) :: dm_aer(ncol, pver, nmodes)  ! per-mode wet number-mode diameter [cm]
      real(kind_phys) :: r_cm, srx_sad(ncol, pver), srx_chk(ncol, pver)

      errmsg = ''
      errflg = 0

      ! Find MAM properties and state from aerosol instances
      aero_props => null()
      aero_state_obj => null()
      do iaermod = 1, aerosol_instances_get_num_models()
        aero_props => aerosol_instances_get_props(iaermod, 0)
        if (associated(aero_props)) then
          if (aero_props%model_is('MAM')) then
            aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
            exit
          end if
        end if
        aero_props => null()
      end do

      if (.not. associated(aero_props) .or. &
          .not. associated(aero_state_obj)) then
         errflg = 1
         errmsg = 'mam_surfarea_diagnostics_run: MAM aerosol model not found in aerosol_instances'
         return
      end if

      ! Tropospheric levels only, as CAM's aero_model_surfarea:
      ! from the level below the tropopause down to the surface.
      beglev(:ncol) = troplev(:ncol) + 1
      endlev(:ncol) = pver

      ! The modal implementation of surf_area_dens ignores relhum (it uses the
      ! wet number-mode diameter from water uptake); the argument exists for
      ! other aerosol representations.
      relhum(:,:) = 0.0_kind_phys

      call aero_state_obj%surf_area_dens(aero_props, sad_chem_spec_types, ncol, pver, &
           beglev, endlev, relhum, pmid, t, pi, sad, reff, sfc, dm_aer)

      call history_out_field('SAD_TROP',  sad)
      call history_out_field('REFF_TROP', reff)
      do m = 1, nmodes
         call history_out_field(trim(sad_name(m)),   sfc(:,:,m))
         call history_out_field(trim(dmaer_name(m)), dm_aer(:,:,m))
      end do

      ! Identity check of each MICM surface reaction's rate parameters against
      ! the MAM modes: 4 pi r^2 N (m2/m3 -> cm2/cm3 is x1e-2) must equal one
      ! mode's surface area density with r (m -> cm is x1e2) half that mode's
      ! wet number-mode diameter; zero where the modes have no surface.
      do j = 1, nsurf
         do k = 1, pver
            do i = 1, ncol
               r_cm = rate_parameters(i,k,surf_iradius(j)) * 1.0e2_kind_phys
               srx_sad(i,k) = 4.0_kind_phys * pi * rate_parameters(i,k,surf_iradius(j))**2 &
                              * rate_parameters(i,k,surf_inumber(j)) * 1.0e-2_kind_phys
               if (all(sfc(i,k,:) == 0.0_kind_phys)) then
                  if (srx_sad(i,k) == 0.0_kind_phys) then
                     srx_chk(i,k) = 0.0_kind_phys
                  else
                     srx_chk(i,k) = -1.0_kind_phys
                  end if
               else
                  srx_chk(i,k) = -1.0_kind_phys
                  do m = 1, nmodes
                     if (sfc(i,k,m) > 0.0_kind_phys .and. &
                         abs(srx_sad(i,k) - sfc(i,k,m)) <= match_rtol * sfc(i,k,m) .and. &
                         abs(r_cm - 0.5_kind_phys * dm_aer(i,k,m)) <= match_rtol * 0.5_kind_phys * dm_aer(i,k,m)) then
                        srx_chk(i,k) = real(m, kind_phys)
                        exit
                     end if
                  end do
               end if
            end do
         end do
         call history_out_field(trim(srx_sad_name(j)), srx_sad)
         call history_out_field(trim(srx_chk_name(j)), srx_chk)
      end do

   end subroutine mam_surfarea_diagnostics_run

end module mam_surfarea_diagnostics
