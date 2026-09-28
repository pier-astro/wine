# Run with Tcl 8.6+: tclsh procedures/test_wine_profile.tcl
source [file join [file dirname [info script]] wine_profile.tcl]

proc _wine_profile_trial {ctx angle} {
    return [expr {($angle - 45.0) * ($angle - 45.0) / 100.0}]
}

set ctx [dict create kind uincl \
    tin [dict create free 0 spec {40 1 0 0 90 90}] \
    tout [dict create free 0 spec {80 1 0 0 90 90}]]
set points {{45 0}}
set lower [_wine_profile_side $ctx 45.0 0.0 1.0 1.0 -1 points]
set upper [_wine_profile_side $ctx 45.0 0.0 1.0 1.0 1 points]
if {$lower ne {boundary 40}} {
    error "expected a boundary-limited lower side, got $lower"
}
if {[lindex $upper 0] ne "ok" ||
    abs([lindex $upper 1] - 55.0) > 0.01} {
    error "expected an upper delta-stat crossing at 55 deg, got $upper"
}
proc _wine_profile_trial {ctx angle} {
    if {$angle == 45.0} {return 0.0}
    return [expr {$angle < 55.0 ? 0.95 : 1.2}]
}
set points {{45 0}}
set approximate [_wine_profile_side $ctx 45.0 0.0 1.0 1.0 1 points]
if {[lindex $approximate 0] ne "approximate" ||
    abs([lindex $approximate 1] - 55.0) > 0.01} {
    error "expected an approximate crossing without warning, got $approximate"
}
set points {{45 0}}
set unstable [_wine_profile_side $ctx 45.0 0.0 1.0 0.1 1 points]
if {[lindex $unstable 0] ne "unstable" ||
    abs([lindex $unstable 1] - 55.0) > 0.001} {
    error "expected a tolerance-exceeding jump near 55 deg, got $unstable"
}
set marker [file tempfile temporary]
close $marker
file delete $temporary
file mkdir $temporary
set ::_wine_profile_failed_trial [file join $temporary failed-trial.xcm]
set ::_wine_profile_output_dir $temporary
proc save {args} {set ::saved $args}
_wine_profile_check_improvement 40.0 9.5 10.0 1.0
if {[info exists ::saved]} {error "insignificant improvement saved a candidate"}
if {![catch {_wine_profile_check_improvement 40.0 8.5 10.0 1.0} message options] ||
    [dict get $options -errorcode] ne \
        [list WINE BETTER_FIT 40.0 8.5 1.5 1.0 better-fit.xcm] ||
    $::saved ne [list all [file join $temporary better-fit.xcm]]} {
    error "significant improvement did not produce an actionable diagnostic: $message"
}
set existing [open [file join $temporary better-fit.xcm] w]
puts $existing original
close $existing
catch {_wine_profile_check_improvement 40.0 8.5 10.0 1.0} message
set input [open [file join $temporary better-fit.xcm] r]
set original [read $input]
close $input
if {$::saved ne [list all [file join $temporary better-fit-2.xcm]] ||
    $original ne "original\n"} {
    error "existing better-fit.xcm was not preserved"
}
set parsed [_wine_profile_arguments 1. 2 toler 0.3 -fit-delta 0.001]
if {$parsed ne {2 1. {} 0.001 0.3}} {error "unexpected profile arguments: $parsed"}
if {![catch {_wine_profile_arguments 1. 2 -toler 0.3}]} {
    error "obsolete -toler spelling was accepted"
}
rename puts _original_puts
proc puts {line} {lappend ::printed $line}
set ::printed {}
_wine_profile_present 43 1.0 43.2073765674 \
    {approximate 36.624501627} {approximate 50.283456789}
if {$::xspec_tclout ne {36.624501627 50.283456789 FFFF} ||
    [string match {*approximate*} [join $::printed]] ||
    ![string match {*36.62*50.28*} [lindex $::printed 2]]} {
    error "approximate profile report lost precision or exposed a quiet status"
}
set ::printed {}
_wine_profile_present 43 1.0 45.0 {boundary 40.0} {unstable 55.123456789}
if {$::xspec_tclout ne {0 55.123456789 TFFT} ||
    ![string match {*0*55.12*} [lindex $::printed 2]] ||
    [llength $::printed] != 5} {
    error "boundary or unstable profile report is incorrect"
}
rename puts {}
rename _original_puts puts
source [file join [file dirname [info script]] wine_inclination.tcl]
proc tclout {option args} {
    upvar 1 xspec_tclout output
    if {$option ne "pinfo"} {error "unexpected tclout $option"}
    set output {incl deg}
}
proc _wine_cov_operand {number args} {return {43.2073765674 1}}
proc _wine_covariance_data {} {return {{1} {0.123456789}}}
rename puts _original_puts
proc puts {line} {lappend ::printed $line}
set ::printed {}
set sigma [wineinclcov 1]
if {$sigma ne $::xspec_tclout || $sigma eq [format %.2f $sigma] ||
    ![string match {*43.21*} [lindex $::printed 0]]} {
    error "covariance return or display precision is incorrect"
}
rename puts {}
rename _original_puts puts
file delete -force $temporary
puts "WINE profile boundary, tolerance, report and diagnostic checks passed"
