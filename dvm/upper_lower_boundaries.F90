#include "fabm_driver.h"

module dvm_upper_lower_boundaries

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries

    type (type_dependency_id)                       :: id_par, id_parmean, id_migrator_food, id_depth 
    type (type_horizontal_dependency_id)            :: id_par0, id_parmean0, id_light_present0, id_nhours, id_migrator_food0, id_topo  
    type (type_horizontal_diagnostic_variable_id)   :: id_nhours_out
    type (type_diagnostic_variable_id)              :: id_present

    contains
        procedure :: initialize
        procedure :: do_surface
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_upper_lower_boundaries), intent(inout), target :: self
        integer, intent(in)                                  :: configunit
        !real(rk) :: par, par0, parmean, parmean0

        call self%register_diagnostic_variable(self%id_present,'migrator_presence','-','migrators are present here')

        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_parmean0,temporal_mean(self%id_par0,period=86400._rk,resolution=3600._rk))
        call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))
        call self%register_dependency(self%id_light_present0,'light_presence','-','light is available at the surface')
        call self%register_dependency(self%id_nhours,temporal_mean(self%id_light_present0,period=86400._rk,resolution=3600._rk))
        call self%register_dependency(self%id_migrator_food,'migrator_food','mgC/m3','food availability for the migrators')
        call self%register_dependency(self%id_migrator_food0,vertical_integral(self%id_migrator_food))
        call self%register_dependency(self%id_depth,standard_variables%pressure)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth )

        call self%register_diagnostic_variable(self%id_nhours_out,'nhours','-','number of daylight hours',source=source_do_surface)
    end subroutine initialize

    subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

        class (type_upper_lower_boundaries),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_SURFACE_

        real(rk) :: nhours
        _HORIZONTAL_LOOP_BEGIN_

            _GET_SURFACE_(self%id_nhours,nhours)
            
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_nhours_out, min(24.0_rk, max(0.0_rk,nhours * 86400_rk)))

        _HORIZONTAL_LOOP_END_

    end subroutine do_surface

    subroutine do(self, _ARGUMENTS_DO_)

        class (type_upper_lower_boundaries), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: par, par0, parmean, parmean0, nhours, food
        real(rk) :: parlog, par0log, parmeanlog, parmean0log
        real(rk) :: depth
        real(rk) :: upper_presence, lower_presence
        real(rk) :: topo
    
        _GET_SURFACE_(self%id_parmean0,parmean0)
        _GET_SURFACE_(self%id_par0,par0)
        _GET_SURFACE_(self%id_nhours,nhours)
        _GET_SURFACE_(self%id_migrator_food0,food)
        _GET_HORIZONTAL_(self%id_topo,topo)

        nhours = min(24.0_rk, max(0.0_rk,nhours * 86400_rk))
        par0log = max(-20.0_rk, log10(par0))
        parmean0log = max(-20.0_rk, log10(parmean0))

        _LOOP_BEGIN_

            _GET_(self%id_par,par)
            _GET_(self%id_parmean,parmean)
            _GET_(self%id_depth,depth)

            parlog = max(-20.0_rk, log10(par))
            parmeanlog = max(-20.0_rk, log10(parmean))

            ! SPECIFY THE POSSIBLE LOCATIONS OF HIGH MIGRATOR CONCENTRATION !

            ! There are 4 cases
            ! 1. Winter Arctic night (surface parmean < 1E-10)
            ! 2. Summer Arctic night (number of daylight hours > 23.9)
            ! 3. Normal day cycle day time
            ! 4. Normal day cycle night time

            ! CASE 1
            upper_presence = 0.0_rk
            lower_presence = 0.0_rk

            if (parmean0 < 1E-10_rk) then
                upper_presence = 0.0_rk
                lower_presence = 0.0_rk
                
                ! Calculate possibilities above the lower boundary
                if (nhours <= 0.30_rk) then
                    if (par0log <= -14.04_rk) then
                        if (food <= 561.66_rk) then
                            if (depth < 164.48_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        else
                            if (depth < 218.94_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        end if
                    else
                        if (par0log <= -8.56_rk) then
                            if (depth < 271.85_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        else
                            if (depth < 398.30_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        end if
                    end if
                else
                    if (food <= 758.34_rk) then
                        if (food <= 758.20_rk) then
                            if (depth < 385.40_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        else
                            if (depth < 470.19_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        end if
                    else
                        if (par0log <= -9.07_rk) then
                            if (depth < 313.93_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        else
                            if (depth < 444.43_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        end if
                    end if
                end if
                
                ! Set diagnostic based on presence
                if (upper_presence + lower_presence > 0.9_rk) then
                    _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                else
                    _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                end if
                
            else

                ! CASE 2
                if (nhours > 23.9_rk) then
                    ! there is an upper and a lower light boundary
                    ! first calculate possibilities above the lower boundary

                    ! Initialize presence variables
                    upper_presence = 0.0_rk
                    lower_presence = 0.0_rk
                    
                    ! Lowerlight Rules
                    if (parmean0log <= 1.22_rk) then
                        if (food <= 785.72_rk) then
                            if (food <= 772.73_rk) then
                                if (parmean0log <= 1.15_rk) then
                                    if (parmeanlog > -15.97_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -18.87_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 0.36_rk) then
                                    if (parmeanlog > -18.20_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -19.27_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (food <= 799.12_rk) then
                                if (par0log <= 0.43_rk) then
                                    if (parmeanlog > -9.16_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -13.40_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (parmean0log <= 1.10_rk) then
                                    if (parmeanlog > -17.41_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -14.79_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= 1.61_rk) then
                            if (food <= 599.12_rk) then
                                if (food <= 570.42_rk) then
                                    if (parmeanlog > -13.88_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -15.96_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (parmean0log <= 1.24_rk) then
                                    if (parmeanlog > -11.76_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -13.75_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= 1.55_rk) then
                                if (par0log <= 1.31_rk) then
                                    if (parmeanlog > -8.54_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -7.68_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 1.73_rk) then
                                    if (parmeanlog > -9.31_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -8.40_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    end if
                    
                    ! Upperlight Rules
                    if (food <= 597.60_rk) then
                        if (par0log <= 0.98_rk) then
                            if (par0log <= 0.45_rk) then
                                if (par0log <= 0.41_rk) then
                                    if (parlog < -5.63_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -5.26_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 567.18_rk) then
                                    if (parlog < -8.08_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.80_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 1.15_rk) then
                                if (par0log <= 1.04_rk) then
                                    if (parlog < -8.43_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -9.78_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 571.23_rk) then
                                    if (parlog < -8.96_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.46_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= 0.92_rk) then
                            if (par0log <= 0.91_rk) then
                                if (food <= 790.53_rk) then
                                    if (parlog < -9.61_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -7.40_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 780.74_rk) then
                                    if (parlog < -10.54_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -9.88_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 1.42_rk) then
                                if (par0log <= 1.04_rk) then
                                    if (parlog < -4.90_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.21_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (parmean0log <= 1.59_rk) then
                                    if (parlog < -3.31_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -2.61_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    end if
                                        
                    ! Set diagnostic based on presence
                    if (upper_presence + lower_presence > 1.0_rk) then
                        _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                    else
                        if (upper_presence > 0.9_rk .and. depth >= max(topo - 20.0_rk, 0.0_rk) ) then 
                            _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                        else 
                            _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                        end if
                    end if
                else

                    ! CASE 3
                    if (par0 > 1E-5_rk) then
                        ! there is an upper and a lower light boundary
                        ! first calculate possibilities above the lower boundary
                        
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Lowerlight Rules
                        if (parmean0log <= 0.86_rk) then
                            if (food <= 562.62_rk) then
                                if (par0log <= 0.07_rk) then
                                    if (parmeanlog > -12.87_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (food <= 562.47_rk) then
                                        if (parmeanlog > -16.12_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -15.37_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (food <= 781.60_rk) then
                                    if (nhours <= 14.59_rk) then
                                        if (parmeanlog > -18.26_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -20.21_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (par0log <= 0.28_rk) then
                                        if (parmeanlog > -13.98_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -17.26_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        else
                            if (nhours <= 21.58_rk) then
                                if (par0log <= 0.56_rk) then
                                    if (par0log <= 0.29_rk) then
                                        if (parmeanlog > -13.67_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -14.41_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (nhours <= 20.45_rk) then
                                        if (parmeanlog > -15.82_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -16.70_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (par0log <= 1.29_rk) then
                                    if (par0log <= 0.79_rk) then
                                        if (parmeanlog > -13.74_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -14.37_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (par0log <= 1.52_rk) then
                                        if (parmeanlog > -14.94_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -15.36_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        end if
                        
                        ! Upperlight Rules
                        if (parmean0log <= 1.07_rk) then
                            if (nhours <= 17.37_rk) then
                                if (food <= 562.45_rk) then
                                    if (parlog < -9.23_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.37_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 762.12_rk) then
                                    if (parlog < -10.16_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.62_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= 1.29_rk) then
                                if (par0log <= 1.01_rk) then
                                    if (parlog < -5.80_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.20_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 1.48_rk) then
                                    if (parlog < -6.67_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.90_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                        
                        ! Set diagnostic based on presence
                        if (upper_presence + lower_presence > 1.0_rk) then
                            _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                        else
                            if (upper_presence > 0.9_rk .and. depth >= max(topo - 20.0_rk, 0.0_rk) ) then 
                                _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                            else 
                                _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                            end if
                        end if
                        
                    else
                        ! CASE 4
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Calculate possibilities above the lower boundary
                        if (food <= 661.97_rk) then
                            if (par0log <= -8.41_rk) then
                                if (food <= 562.06_rk) then
                                    if (food <= 562.06_rk) then
                                        if (parmeanlog > -9.16_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -8.40_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (food <= 562.08_rk) then
                                        if (parmeanlog > -11.89_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -10.26_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (nhours <= 20.52_rk) then
                                    if (nhours <= 20.09_rk) then
                                        if (parmeanlog > -14.09_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -12.20_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (food <= 563.64_rk) then
                                        if (parmeanlog > -15.38_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -16.43_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 0.60_rk) then
                                if (par0log <= -1.07_rk) then
                                    if (food <= 760.25_rk) then
                                        if (parmeanlog > -19.94_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -15.65_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (par0log <= -0.29_rk) then
                                        if (parmeanlog > -20.17_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -19.04_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (nhours <= 16.94_rk) then
                                    if (parmeanlog > -9.78_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (nhours <= 20.52_rk) then
                                        if (parmeanlog > -13.61_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -11.00_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        end if
                        
                        ! Set diagnostic based on presence
                        if (upper_presence + lower_presence > 0.9_rk) then
                            _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                        else
                            _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                        end if
                        
                    end if
                end if

            end if

            ! This should ensure that each point at least receives the 0.0_rk value
            if (upper_presence + lower_presence < 0.9_rk) then
                _SET_DIAGNOSTIC_(self%id_present,0.0_rk)     
            end if 
            ! 


            ! ------------------------------------------------------------- !


        _LOOP_END_

    end subroutine do

end module
