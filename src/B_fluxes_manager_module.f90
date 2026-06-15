module fluxes_manager_module
use mesh_module, only: mesh_t
use field_manager_module, only: field_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use solver_config_module, only: db_solver_config_t
use fluxes_calculation_module
use ADdiff_fluxes_calculation_module
use stability_operator_module, only: get_J_mtrx
implicit none
!=======================================================================
!==================== FLUXES MANAGER DATA STRUCTURE ====================
!=======================================================================
type, public :: flux_manager_t
    type(mesh_t), pointer :: mesh => null()
    type(field_manager_t), pointer :: fldsm => null()
    type(phys_prop_manager_t), pointer :: ppm => null()
    type(db_solver_config_t), pointer :: settings => null()

contains
    procedure :: initialize => flxsm_initialize
    procedure :: compute_fluxes_cmprs => flxsm_compute_fluxes_cmprs
    procedure :: compute_ADdiff_fluxes_cmprs => flxsm_compute_ADdiff_fluxes_cmprs
	procedure :: cons_to_prim => flxsm_cons_to_prim_cmprs
    procedure :: compute_inviscid => flxsm_compute_inviscid_cmprs
    procedure :: ADdiff_compute_inviscid => flxsm_compute_ADdiff_inviscid_cmprs
    procedure :: compute_viscous  => flxsm_compute_viscous_cmprs
    procedure :: ADdiff_compute_viscous  => flxsm_compute_ADdiff_viscous_cmprs
end type

contains
!=======================================================================
!======================= FLUXES MANAGER METHODS ========================
!=======================================================================
	subroutine flxsm_initialize(this, mesh, fldsm, ppm, settings)
		implicit none
		class(flux_manager_t), intent(inout) :: this
		type(mesh_t), intent(in), target:: mesh
		type(field_manager_t), intent(in), target :: fldsm
		type(phys_prop_manager_t), intent(in), target :: ppm
		type(db_solver_config_t), intent(in), target :: settings
		
		this%mesh => mesh
		this%fldsm => fldsm
		this%ppm => ppm
		this%settings => settings		
	end subroutine

	subroutine flxsm_compute_fluxes_cmprs(this, vars_idx, fluxes_idx)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, fluxes_idx
		
		!=== INVISCID FLUXES: ===
		call this%compute_inviscid(vars_idx, fluxes_idx, .true.)
		
		!=== VISCOUS FLUXES: ===
		if (this%settings%MODEL .NE. 1) then
			call this%compute_viscous(vars_idx, fluxes_idx, .false.)
		end if
		
		!=== CONSERVATIVE TO PRIMITIVE: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			call this%cons_to_prim(vars_idx, fluxes_idx)
		end if
	end subroutine
	
	subroutine flxsm_compute_ADdiff_fluxes_cmprs(this, vars_idx, dvars_idx, dfluxes_idx)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, dvars_idx, dfluxes_idx
		
		!=== INVISCID FLUXES: ===
		call this%ADdiff_compute_inviscid(vars_idx, dvars_idx, dfluxes_idx, .true.)
		
		!=== VISCOUS FLUXES: ===
		if (this%settings%MODEL .NE. 1) then
			call this%ADdiff_compute_viscous(vars_idx, dvars_idx, dfluxes_idx, .false.)
		end if
		
		!=== CONSERVATIVE TO PRIMITIVE: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			call this%cons_to_prim(vars_idx, dfluxes_idx)
		end if
	end subroutine
	
	
	subroutine flxsm_compute_inviscid_cmprs(this, vars_idx, fluxes_idx, is_vanish)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, fluxes_idx
		logical, intent(in) :: is_vanish
		
		real(8), pointer, contiguous :: Q(:, :), gradQ(:, :), Fluxes(:, :)
		
		Q => this%fldsm%registry(vars_idx)%values
		gradQ => this%fldsm%registry(vars_idx)%grad
		Fluxes => this%fldsm%registry(fluxes_idx)%values
		
		if (is_vanish) Fluxes = 0.d0
										   
		call compute_inviscid_fluxes_cmprs(this%mesh%dim,&
										   this%mesh%nfaces,&
										   this%mesh%nbfaces,&
									       this%mesh%face_left_cell,&
									       this%mesh%face_right_cell,&
										   this%mesh%face_normal,&
										   this%mesh%face_area,&
										   this%mesh%face_center,&
										   this%mesh%cell_center,&
										   Q, gradQ, Fluxes,&
										   this%ppm%k, this%ppm%R_gas, this%ppm%cp, this%ppm%cv,&
										   this%settings%ORDER, this%settings%LIMITER,&
										   this%settings%REIMAN_SOLVER_TYPE,&
										   this%settings%USE_GHOST_CELLS)
										   								   
	end subroutine
	
	subroutine flxsm_compute_ADdiff_inviscid_cmprs(this, vars_idx, dvars_idx, dfluxes_idx, is_vanish)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, dvars_idx, dfluxes_idx
		logical, intent(in) :: is_vanish
		
		real(8), pointer, contiguous :: Q(:, :), dQ(:, :), gradQ(:, :),&
		                                dgradQ(:, :), dFluxes(:, :)
		
		Q => this%fldsm%registry(vars_idx)%values
		gradQ => this%fldsm%registry(vars_idx)%grad
		dQ => this%fldsm%registry(dvars_idx)%values
		dgradQ => this%fldsm%registry(dvars_idx)%grad
		dFluxes => this%fldsm%registry(dfluxes_idx)%values
		
		if (is_vanish) dFluxes = 0.d0
		call COMPUTE_INVISCID_FLUXES_CMPRS_D(this%mesh%dim,&
										     this%mesh%nfaces,&
										     this%mesh%nbfaces,&
										     this%mesh%face_left_cell,&
									         this%mesh%face_right_cell,&
									         this%mesh%face_normal,&
										     this%mesh%face_area,&
										     this%mesh%face_center,&
										     this%mesh%cell_center,&
											 Q, dQ, gradQ, dgradQ, dFluxes,&
											 this%ppm%k, this%ppm%R_gas, this%ppm%cp, this%ppm%cv,&
										     this%settings%ORDER, this%settings%LIMITER,&
										     this%settings%REIMAN_SOLVER_TYPE,&
										     this%settings%USE_GHOST_CELLS)										   								   
	end subroutine
	
	
	subroutine flxsm_compute_viscous_cmprs(this, vars_idx, fluxes_idx, is_vanish)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, fluxes_idx
		logical, intent(in) :: is_vanish
		
		real(8), pointer, contiguous :: Q(:, :), gradQ(:, :), Fluxes(:, :)
		
		Q => this%fldsm%registry(vars_idx)%values
		gradQ => this%fldsm%registry(vars_idx)%grad
		Fluxes => this%fldsm%registry(fluxes_idx)%values
		
		if (is_vanish) Fluxes = 0.d0												 
												 
		call compute_viscous_fluxes_cmprs(this%mesh%dim,&
										  this%mesh%nfaces,&
										  this%mesh%nbfaces,&
										  this%mesh%face_left_cell,&
										  this%mesh%face_right_cell,&
										  this%mesh%face_normal,&
										  this%mesh%face_area,&
										  this%mesh%face_weight,&
										  this%mesh%face_center,&
										  this%mesh%cell_center,&
										  Q, gradQ, Fluxes,&
										  this%ppm%k, this%ppm%R_gas, this%ppm%cp,&
										  this%ppm%cv, this%ppm%Pr,&
										  this%settings%USE_GHOST_CELLS) 
	end subroutine
	
	subroutine flxsm_compute_ADdiff_viscous_cmprs(this, vars_idx, dvars_idx, dfluxes_idx, is_vanish)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, dvars_idx, dfluxes_idx
		logical, intent(in) :: is_vanish
		
		real(8), pointer, contiguous :: Q(:, :), dQ(:, :), gradQ(:, :),&
		                                dgradQ(:, :), dFluxes(:, :)
		
		Q => this%fldsm%registry(vars_idx)%values
		gradQ => this%fldsm%registry(vars_idx)%grad
		dQ => this%fldsm%registry(dvars_idx)%values
		dgradQ => this%fldsm%registry(dvars_idx)%grad
		dFluxes => this%fldsm%registry(dfluxes_idx)%values
		
		if (is_vanish) dFluxes = 0.d0
		
		call COMPUTE_VISCOUS_FLUXES_CMPRS_D(this%mesh%dim,&
										    this%mesh%nfaces,&
										    this%mesh%nbfaces,&
										    this%mesh%face_left_cell,&
									        this%mesh%face_right_cell,&
									        this%mesh%face_normal,&
										    this%mesh%face_area,&
										    this%mesh%face_weight,&
										    this%mesh%face_center,&
										    this%mesh%cell_center,&
										    Q, dQ, gradQ, dgradQ, dFluxes,&
										    this%ppm%k, this%ppm%R_gas, this%ppm%cp,&
										    this%ppm%cv, this%ppm%Pr,&
										    this%settings%USE_GHOST_CELLS)									   								   
	end subroutine
	
	
	subroutine flxsm_cons_to_prim_cmprs(this, vars_idx, fluxes_idx)
		implicit none
		class(flux_manager_t), intent(inout), target :: this
		integer, intent(in) :: vars_idx, fluxes_idx
		
		real(8), pointer, contiguous :: Q(:, :), Fluxes(:, :)
		
		Q => this%fldsm%registry(vars_idx)%values
		Fluxes => this%fldsm%registry(fluxes_idx)%values
		
		call cons_fluxes_to_prim_fluxes(this%mesh%dim,&
										this%mesh%ncells,&
										Q, Fluxes,&
										this%ppm%k, this%ppm%R_gas,&
										this%ppm%cv, this%ppm%cp)
	end subroutine

	
	

	SUBROUTINE CONS_TO_PRIM_D(dim, ncells, cv, r_gas, u, ud, q, qd)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: dim, ncells
		REAL*8, INTENT(IN) :: cv, r_gas
		REAL*8, INTENT(IN), CONTIGUOUS :: u(:, :)
		REAL*8, INTENT(IN), CONTIGUOUS :: ud(:, :)
		REAL*8, INTENT(IN), CONTIGUOUS :: q(:, :)
		REAL*8, INTENT(INOUT), CONTIGUOUS :: qd(:, :)
		REAL*8 :: v2d
		INTEGER :: cell_idx, d, v
		REAL*8 :: temp0
	!recomputing primitive variables:
		DO cell_idx=1,ncells
	!velocity:
		  DO d=2,dim+1
			temp0 = u(d, cell_idx)/u(1, cell_idx)
			qd(d, cell_idx) = (ud(d, cell_idx)-temp0*ud(1, cell_idx))/u(1, &
	&         cell_idx)
		  END DO
	!temperature:
		  v2d = 0.0_8
		  DO v=2,dim+1
			v2d = v2d + 2*q(v, cell_idx)*qd(v, cell_idx)
		  END DO
		  temp0 = u(dim+2, cell_idx)/u(1, cell_idx)
		  qd(dim+2, cell_idx) = ((ud(dim+2, cell_idx)-temp0*ud(1, cell_idx))&
	&       /u(1, cell_idx)-0.5d0*v2d)/cv
	!pressure:
		  qd(1, cell_idx) = r_gas*(q(dim+2, cell_idx)*ud(1, cell_idx)+u(1, &
	&       cell_idx)*qd(dim+2, cell_idx))
		END DO
	  END SUBROUTINE CONS_TO_PRIM_D
end module
