#include "fabm_driver.h"

module dvm_upper_lower_boundaries

use fabm_types
use fabm_expressions

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries

    type (type_dependency_id)                       :: id_par, id_parmean, id_migrator_food, id_depth 
    type (type_horizontal_dependency_id)            :: id_par0, id_parmean0, id_light_present0, id_nhours, id_migrator_food0  
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
    
        _GET_SURFACE_(self%id_parmean0,parmean0)
        _GET_SURFACE_(self%id_par0,par0)
        _GET_SURFACE_(self%id_nhours,nhours)
        _GET_SURFACE_(self%id_migrator_food0,food)

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
            upper_presence = 0.0_rk
            lower_presence = 0.0_rk

            ! CASE 1
            if (parmean0 < 1E-10_rk) then
                upper_presence = 0.0_rk
                lower_presence = 0.0_rk
                if (nhours <= 0.3_rk) then
                    if (par0log <= -14.04_rk) then
                        if (food <= 757.41_rk) then
                            if (depth < 198.67_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        else
                            if (depth < 236.16_rk) then
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
                    if (food <= 995.4_rk) then
                        if (nhours <= 7.15_rk) then
                            if (depth < 446.33_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if                            
                        else
                            if (depth < 344.38_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        end if
                    else
                        if (food <= 997.81_rk) then
                            if (depth < 256.12_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        else
                            if (depth < 172.79_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                        end if
                    end if                                                          

                end if
                if (upper_presence + lower_presence > 0.9_rk) then ! a band in the water column should be 1.0
                    _SET_DIAGNOSTIC_(self%id_present,1.0_rk)       ! lower presence is ignored in this case 
                else                                               ! since migrator will search the surface by default
                    _SET_DIAGNOSTIC_(self%id_present,0.0_rk)
                end if 
            else

                ! CASE 2
                if (nhours > 23.9_rk) then
                    ! there is an upper and a lower light boundary
                    ! first calculate possibilities above the lower boundary
                    upper_presence = 0.0_rk
                    lower_presence = 0.0_rk
                    if (parmean0log <= 1.22_rk) then
                        if (food <= 781.76_rk) then
                            if (food <= 163.16_rk) then
                                if (parlog > -15.76_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            else
                                if (parlog > -19.05_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            end if
                        else
                            if (food <= 791.07_rk) then
                                if (parlog > -13.17_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            else
                                if (parlog > -16.92_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= 1.61_rk) then
                            if (parmean0log <= 1.24_rk) then
                                if (parlog > -11.92_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            else
                                if (parlog > -14.80_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            end if
                        else
                            if (food <= 2109.83_rk) then
                                if (parlog > -8.36_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            else
                                if (parlog > -8.90_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            end if
                        end if
                    end if  

                    ! now calculate possibilities below the upper boundary
                    if (food <= 1298.09_rk) then
                        if (parmean0log <= 1.49_rk) then
                            if (par0log <= 0.91_rk) then
                                if (parlog < -8.23_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            else
                                if (parlog < -9.41_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            end if
                        else
                            if (parmean0log <= 1.55_rk) then
                                if (parlog < -3.27_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            else
                                if (parlog < -8.42_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= 1.48_rk) then
                            if (food <= 2076.14_rk) then
                                if (parlog < -4.25_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            else
                                if (parlog < -6.18_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            end if
                        else
                            if (food <= 2119.90_rk) then
                                if (parlog < -2.24_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            else
                                if (parlog < -2.60_rk) then
                                    lower_presence = 1.0_rk
                                else
                                    lower_presence = 0.0_rk
                                end if
                            end if
                        end if
                    end if
                    
                    if (upper_presence + lower_presence > 1.0_rk) then ! a band in the water column should be 2.0
                        _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                    else
                        _SET_DIAGNOSTIC_(self%id_present,0.0_rk)
                    end if 
                    
                else

                    ! CASE 3
                    if (par0 > 1E-5_rk) then
                        ! there is an upper and a lower light boundary
                        ! first calculate possibilities above the lower boundary
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        if (food <= 458.55_rk) then
                            if (parmean0log <= 0.39_rk) then
                                if (par0log <= 0.07_rk) then
                                    if (parlog > -12.53_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog > -15.40_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 154.52_rk) then
                                    if (parlog > -17.91_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog > -15.84_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 0.60_rk) then
                                if (par0log <= 0.61_rk) then
                                    if (parlog > -19.54_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog > -19.98_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 2098.06_rk) then
                                    if (parlog > -17.21_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog > -14.92_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if  

                        ! now calculate possibilities below the upper boundary            
                        if (parmean0log <= 1.07_rk) then
                            if (nhours <= 17.37_rk) then
                                if (par0log <= 0.30_rk) then
                                    if (parlog < -7.61_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.76_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 762.30_rk) then
                                    if (parlog < -10.18_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.75_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= 1.29_rk) then
                                if (par0log <= 1.01_rk) then
                                    if (parlog < -5.93_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.35_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 1.48_rk) then
                                    if (parlog < -6.82_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -7.05_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if

                        if (upper_presence + lower_presence > 1.0_rk) then ! a band in the water column should be 2.0
                            _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                        else
                            _SET_DIAGNOSTIC_(self%id_present,0.0_rk)
                        end if

                    else
                        ! CASE 4
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        if (food <= 155.33_rk) then
                            if (par0log <= -8.41_rk) then
                                if (food <= 151.89_rk) then
                                    if (parmeanlog > -13.07_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -9.9_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (nhours <= 8.93_rk) then
                                    if (parmeanlog > -11.93_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -14.44_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 0.6_rk) then
                                if (par0log <= -1.07_rk) then
                                    if (parmeanlog > -17.89_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog >= -20.0_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 155.71_rk) then
                                    if (parmeanlog >= -15.35_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog >= -11.71_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if                                                                                               
                        end if
                        if (upper_presence + lower_presence > 0.9_rk) then ! a band in the water column should be 1.0
                            _SET_DIAGNOSTIC_(self%id_present,1.0_rk)       ! lower presence is ignored in this case 
                        else                                               ! since migrator will search the surface by default
                            _SET_DIAGNOSTIC_(self%id_present,0.0_rk)
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