# Parameter-layout regression for native uniform-width XSPEC geometry.
source [file join [file dirname [info script]] .. procedures wind_geometry.tcl]
rename _wind_out _wind_real_out
proc _wind_out {command args} {
    if {$command eq "compinfo"} {return [list $::model 1 [expr {$::offset+3}]]}
    set index [lindex $args 0]
    set relative [expr {$index-$::offset-1}]
    if {$command eq "pinfo"} {return [list [lindex $::keys $relative] {}]}
    if {$command eq "plink"} {return F}
    if {$command eq "pfree"} {return T}
    if {$command eq "param"} {return [list [lindex $::values $relative] .01 0 0 90 90]}
    error "Unexpected command $command"
}
foreach {model offset keys values} {
    wemlos_Cf 1 {Cf skew incl} {.7 .5 45}
    wuconvlos_Cf 2 {Cf skew incl} {.7 .5 45}
    wuconv 2 {tin tout incl} {20 70 45}
} {
    set context [_wind_context 1 [_wind_options]]
    set actual [_wind_angles $context]
    if {[dict get $context [lindex $keys 0] id] != $offset+1} {error "Wrong geometry offset"}
    if {abs([dict get $actual incl]-45)>1e-10} {error "Wrong inclination"}
    if {[dict get $actual tin]>[dict get $actual incl] || [dict get $actual tout]<[dict get $actual incl]} {
        error "Sightline outside expected cone"
    }
    if {$model ne "wuconv" && abs([dict get $actual Cf]-.7)>1e-12} {error "Wrong covering"}
    puts "$model geometry layout passed"
}
