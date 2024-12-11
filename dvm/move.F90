#include "fabm_driver.h"

module dvm_move

use fabm_types
use fabm_expressions

implicit none 

private 

type, extends(type_base_model), public :: type_move
    type (type_dependency_id)          :: id_random_weights
    type (type_dependency_id)          :: id_thickness
    type (type_state_variable_id)      :: id_target
    type (type_bottom_dependency_id)   :: id_integral
    type (type_bottom_dependency_id)   :: id_integral_random_weights
    type (type_diagnostic_variable_id) :: id_distributed
    real(rk) :: tstep, m, ratioMig
    contains
        procedure ::initialize
        procedure :: do
end type

contains

    subroutine initialize(self, configunit)

        class (type_move), intent(inout), target :: self
        integer,         intent(in)            :: configunit

        call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
        call self%register_dependency(self%id_random_weights,'migrator_random_weights','-','migrators distribution random weights')
        call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')
        call self%register_dependency(self%id_integral_random_weights,'migrator_integral_random_weights','','migrators distribution integral random weights')
        call self%register_diagnostic_variable(self%id_distributed,'migrator_distributed_weights','-','migrators distribution random weights', missing_value=0.0_rk, source=source_do)
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%get_parameter( self%tstep,  'tstep',    'sec',      'time-step in seconds', default=600.0_rk)
        call self%get_parameter( self%ratioMig,  'ratioMig',    '-',      'ratio of moving biomass', default=0.5_rk)

    end subroutine initialize

    subroutine do(self, _ARGUMENTS_DO_)
        class (type_move), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: integral, integral_random_weights, random_weights
        real(rk) :: local, distributed, thickness
        real(rk) :: mortality_switch, local_loss
    
        _GET_BOTTOM_(self%id_integral,integral)
        _GET_BOTTOM_(self%id_integral_random_weights,integral_random_weights)
        !write(*,*)integral
        _LOOP_BEGIN_
           _GET_(self%id_target,local)
           _GET_(self%id_random_weights,random_weights)
           _GET_(self%id_thickness,thickness)
    
           distributed = random_weights / integral_random_weights ! this ensures that total weight distribution is 1, so mass conserved
           !write(*,*)distributed
           _SET_DIAGNOSTIC_(self%id_distributed,distributed)
           !write(*,*)par,mortality_switch,local_loss,local
           _ADD_SOURCE_(self%id_target, (distributed * integral/thickness - local ) * self%ratioMig / self%tstep) 
           !write(*,*)self%dt
    
        _LOOP_END_
     end subroutine do

end module
