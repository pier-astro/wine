# Pure geometry derivatives and covariance, without XSPEC.
source [file join [file dirname [info script]] wind_geometry.tcl]
proc _wind_out {option id} {return [list [dict get $::values $id]]}

foreach kind {generic los off} {
    set shape [expr {$kind eq "los" ? "skew" : "hole"}]
    set view [expr {$kind eq "off" ? "view" : "incl"}]
    set ::values [dict create 1 0.4 2 0.35 3 [expr {$kind eq "off" ? 0.4 : 35.0}]]
    foreach region {0 1} {
        set ctx [dict create kind $kind region $region]
        set coordinates {}
        foreach key [list Cf $shape $view] id {1 2 3} {
            dict set ctx $key [dict create id $id free 1]
            dict set coordinates $key [dict get $::values $id]
        }
        set geometry [_wind_angles $ctx $coordinates]
        foreach quantity {Cf tin tout incl} {
            set gradient [_wind_gradient $ctx $geometry $quantity]
            foreach key [list Cf $shape $view] id {1 2 3} {
                set value [dict get $coordinates $key]
                set step 1e-6
                set plus $coordinates
                set minus $coordinates
                dict set plus $key [expr {$value+$step}]
                dict set minus $key [expr {$value-$step}]
                set numerical [expr {([dict get [_wind_angles $ctx $plus] $quantity]-
                    [dict get [_wind_angles $ctx $minus] $quantity])/(2.0*$step)}]
                if {abs($numerical-[dict get $gradient $id])>1e-5*max(1.0,abs($numerical))} {
                    _wind_fail "Incorrect $kind/$region $quantity derivative with respect to $key"
                }
            }
        }
    }
}
set propagated [_wind_covariance_error [dict create 1 2.0 2 -1.0] {{1 2} {4 1 9}}]
if {abs($propagated-sqrt(21.0))>1e-12} {_wind_fail "Missing covariance cross term"}
puts "WINE analytic geometry derivatives and full-covariance propagation passed"

source [file join [file dirname [info script]] wine_profile.tcl]
proc _wind_test_trial {ctx angle} {return [expr {$angle*$angle}]}
set points {}
set profile [dict create profile_bounds {0.00001 89.99999}]
set lower [_wine_profile_side $profile 0.0 0.0 1.0 1.0 -1 points _wind_test_trial]
set upper [_wine_profile_side $profile 90.0 0.0 1.0 1.0 1 points _wind_test_trial]
if {[lindex $lower 0] ne "boundary" || [lindex $lower 1]!=0 ||
    [lindex $upper 0] ne "boundary" || [lindex $upper 1]!=90} {_wind_fail "Incorrect pegged-center profile boundary"}
puts "WINE shared profile callback and boundary-center checks passed"
