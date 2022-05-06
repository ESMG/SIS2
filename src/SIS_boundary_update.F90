! This file is part of SIS2. See LICENSE.md for the license.

!> Controls where open boundary conditions are applied
module SIS_boundary_update

use ice_grid,                 only : ice_grid_type
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,        only : time_type
use MOM_error_handler,        only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,          only : get_param, log_version, param_file_type, log_param
use MOM_grid,                 only : ocean_grid_type
use MOM_dyn_horgrid,          only : dyn_horgrid_type
use MOM_unit_scaling,         only : unit_scale_type
use SIS_open_boundary,        only : ice_OBC_registry_type, file_ice_OBC_CS
use SIS_open_boundary,        only : register_file_ice_OBC, file_ice_OBC_end
use SIS_open_boundary,        only : ice_OBC_type, update_ice_OBC_segment_data
use SIS_tracer_registry,      only : SIS_tracer_registry_type

implicit none ; private

#include <MOM_memory.h>

public call_ice_OBC_register, ice_OBC_register_end
public update_ice_OBC_data

!> The control structure for the MOM_boundary_update module
type, public :: update_ice_OBC_CS ; private
  logical :: use_files = .false.        !< If true, use external files for the open boundary.
  !>@{ Pointers to the control structures for named OBC specifications
  type(file_ice_OBC_CS), pointer :: file_ice_OBC_CSp => NULL()
  !>@}
end type update_ice_OBC_CS

integer :: id_clock_pass !< A CPU time clock ID

! character(len=40)  :: mdl = "MOM_boundary_update" ! This module's
! name.

contains

!> The following subroutines and associated definitions provide the
!! machinery to register and call the subroutines that initialize
!! open boundary conditions.
subroutine call_ice_OBC_register(param_file, CS, US, OBC, tr_Reg)
  type(param_file_type),     intent(in)   :: param_file !< Parameter file to parse
  type(update_ice_OBC_CS),   pointer      :: CS         !< Control structure for OBCs
  type(unit_scale_type),     intent(in)   :: US         !< A dimensional unit scaling type
  type(ice_OBC_type),        pointer      :: OBC        !< Open boundary structure
  type(SIS_tracer_registry_type), pointer :: tr_Reg     !< Tracer registry.

  ! Local variables
  character(len=200) :: config
  character(len=40)  :: mdl = "SIS_boundary_update" ! This module's name.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  if (associated(CS)) then
    call MOM_error(WARNING, "call_ice_OBC_register called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif


  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "USE_FILE_OBC", CS%use_files, &
                 "If true, use external files for the open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "OBC_USER_CONFIG", config, &
               "A string that sets how the user code is invoked to set open boundary data: \n"//&
               "   USER - user specified", default="none", do_not_log=.true.)

  if (CS%use_files) CS%use_files = &
    register_file_ice_OBC(param_file, CS%file_ice_OBC_CSp, US, &
               OBC%OBC_Reg)

end subroutine call_ice_OBC_register

!> Calls appropriate routine to update the open boundary conditions.
subroutine update_ice_OBC_data(OBC, G, IG, US, CS, Time)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean grid structure
  type(ice_grid_type),                       intent(in)    :: IG   !< Ice vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
! type(thermo_var_ptrs),                     intent(in)    :: tv   !< Thermodynamics structure
! real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h    !< layer thicknesses [H ~> m or kg m-2]
  type(ice_OBC_type),                        pointer       :: OBC  !< Open boundary structure
  type(update_ice_OBC_CS),                   pointer       :: CS   !< Control structure for ice OBCs
  type(time_type),                           intent(in)    :: Time !< Model time

! Something here... with CS%file_ice_OBC_CSp?
! if (CS%use_files) &
!     call update_ice_OBC_segment_data(G, GV, OBC, tv, h, Time)
  if (OBC%needs_IO_for_data) call update_ice_OBC_segment_data(G, IG, US, OBC, Time)

end subroutine update_ice_OBC_data

!> Clean up the OBC registry.
subroutine ice_OBC_register_end(CS)
  type(update_ice_OBC_CS),       pointer    :: CS !< Control structure for OBCs

  if (CS%use_files) call file_ice_OBC_end(CS%file_ice_OBC_CSp)
  if (associated(CS)) deallocate(CS)
end subroutine ice_OBC_register_end

!> \namespace SIS_boundary_update
!! This module updates the open boundary arrays when time-varying.
!! It caused a circular dependency with the tidal_bay and other setups when in
!! SIS_open_boundary.
!!
!! A small fragment of the grid is shown below:
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, tauy
!!    j    x ^ x ^ x   At >:  u, taux
!!    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!!
!! The boundaries always run through q grid points (x).

end module SIS_boundary_update
