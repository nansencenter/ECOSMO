#include "fabm_driver.h"

module dvm_weight_distribution

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_weight_distribution
    type (type_bottom_diagnostic_variable_id) :: id_integral
    type (type_bottom_diagnostic_variable_id) :: id_integral_random_weights
    type (type_diagnostic_variable_id)        :: id_random_weights
    type (type_state_variable_id)             :: id_target
    type (type_dependency_id)                 :: id_present
    type (type_dependency_id)                 :: id_migrator_food
    type (type_dependency_id)                 :: id_thickness
    type (type_horizontal_dependency_id)      :: id_par0 !, id_parmean0

    contains
        procedure :: initialize
        procedure :: do_column
end type

contains

    subroutine initialize(self, configunit)
        class (type_weight_distribution), intent(inout), target :: self
        integer,                     intent(in)            :: configunit
        call self%register_diagnostic_variable(self%id_integral_random_weights,'migrator_integral_random_weights','-','migrators distribution integral random weights', missing_value=0.0_rk, &
            act_as_state_variable=.true., source=source_do_column)
        call self%register_diagnostic_variable(self%id_random_weights,'migrator_random_weights','-','migrators distribution random weights', missing_value=0.0_rk, &
            act_as_state_variable=.true., source=source_do_column)
        call self%register_state_dependency(self%id_target, 'target', '', 'variable to depth-integrate')   
        call self%register_dependency(self%id_present,'present', '-', 'migrators are present here')
        call self%register_dependency(self%id_migrator_food,'migrator_food', 'mgC/m3', 'food availability for the migrators')
        call self%register_diagnostic_variable(self%id_integral, 'integral', '','integral', missing_value=0.0_rk, &
            act_as_state_variable=.true., source=source_do_column)
        call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
!    call self%register_dependency(self%id_parmean0,temporal_mean(self%id_par0,period=86400._rk,resolution=3600._rk))

        call random_seed()
    end subroutine initialize

    subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
        class (type_weight_distribution), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_COLUMN_

        real(rk) :: local, weight, thickness, integral, integral_random, present, local_random, minimum_random, minimum_value
        real(rk) :: totaldepth
        real(rk) :: par0
        real(rk) :: food

        integral = 0.0_rk
        integral_random = 0.0_rk
        totaldepth = 0.0_rk

        _GET_SURFACE_(self%id_par0,par0)

        _VERTICAL_LOOP_BEGIN_

            _GET_(self%id_target,local)
            _GET_(self%id_present,present)
            _GET_(self%id_thickness,thickness)
            _GET_(self%id_migrator_food,food)
            integral = integral + local*thickness
            totaldepth = totaldepth + thickness

            call random_number(local_random)
            call random_number(minimum_random)
       !local_random = 0.1_rk + 0.9_rk * local_random * present
       !local_random = local_random * present
            minimum_value = 0.002 + (0.01 - 0.002) * minimum_random
       !write(*,*)parmean0
    !    if (parmean0 < 1E-4_rk) then
    !         if (totaldepth <= 150.0_rk + minimum_value * 1000.0_rk - 40.0_rk) then ! adds around 50m randomness around 150 meters
    !             local_random = thickness * (minimum_value + (1.0_rk - minimum_value) * local_random)
    !         else
    !             local_random = thickness * minimum_value
    !         end if 
    !    else
            if (par0 <= 1E-18_rk) then 
                local_random = thickness * (minimum_value + (1.0_rk - minimum_value) * local_random * present * food)
            else
                local_random = thickness * (minimum_value + (1.0_rk - minimum_value) * local_random * present)
            end if
       !write(*,*)local_random
            integral_random = integral_random + local_random

            _SET_DIAGNOSTIC_(self%id_random_weights,local_random)

        _VERTICAL_LOOP_END_
        _SET_BOTTOM_DIAGNOSTIC_(self%id_integral,integral)
        _SET_BOTTOM_DIAGNOSTIC_(self%id_integral_random_weights,integral_random)

    end subroutine do_column

end module
