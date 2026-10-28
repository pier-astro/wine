# Physical geometry reporting and isolated angular scans for current WINE.
# Public selectors are component numbers, with optional model NAME/group N.
set ::_wind_source [file normalize [info script]]
if {[file type $::_wind_source] eq "link"} {
    set ::_wind_source [file normalize [file join [file dirname $::_wind_source] \
        [file readlink $::_wind_source]]]
}
set ::_wind_library [file normalize [file join [file dirname $::_wind_source] .. xspec]]

proc _wind_fail {message} {return -code error $message}

proc _wind_out {args} {
    tclout {*}$args
    return $xspec_tclout
}

proc _wind_options {args} {
    set options [dict create model {} group 1]
    if {[llength $args] % 2} {_wind_fail "Expected model NAME and/or group N"}
    foreach {key value} $args {
        if {![dict exists $options $key]} {_wind_fail "Unknown option: $key"}
        dict set options $key $value
    }
    set group [dict get $options group]
    if {![string is integer -strict $group] || $group < 1} {_wind_fail "group must be a positive integer"}
    return $options
}

proc _wind_parameter {id expected {seen {}}} {
    if {$id in $seen} {_wind_fail "Cyclic link at $id"}
    if {[lindex [_wind_out pinfo $id] 0] ne $expected} {_wind_fail "Expected $expected at $id"}
    set link [_wind_out plink $id]
    if {[lindex $link 0] eq "T"} {
        if {![regexp {^T\s+=\s+((?:[A-Za-z_][A-Za-z_0-9]*:)?p?[0-9]+)$} $link -> target]} {
            _wind_fail "Unsupported non-identity geometry link at $id: $link"
        }
        regsub {(^|:)p} $target {\1} target
        return [_wind_parameter $target $expected [linsert $seen end $id]]
    }
    set spec [_wind_out param $id]
    return [dict create id $id spec $spec value [lindex $spec 0] \
        free [expr {[_wind_out pfree $id] eq "T"}]]
}

proc _wind_context {component options} {
    if {![string is integer -strict $component] || $component < 1} {_wind_fail "Expected a positive component number"}
    set prefix [dict get $options model]
    if {$prefix ne ""} {append prefix :}
    lassign [_wind_out compinfo ${prefix}$component [dict get $options group]] name first count
    if {$name ni {windem windemlos windemoff windline windlinelos windlineoff}} {
        _wind_fail "Component $component is not a WINE emission component: $name"
    }
    set kind generic
    if {[string match *los $name]} {set kind los}
    if {[string match *off $name]} {set kind off}
    set shape [expr {$kind eq "los" ? "skew" : "hole"}]
    set view [expr {$kind eq "off" ? "view" : "incl"}]
    set ctx [dict create component $component options $options kind $kind name $name]
    foreach key [list Cf $shape $view] offset {1 2 3} {
        dict set ctx $key [_wind_parameter ${prefix}[expr {$first+$offset}] $key]
    }
    if {$kind eq "off"} {
        set region [_wind_parameter ${prefix}[expr {$first+4}] region]
        if {[dict get $region free] || [dict get $region value] ni {0 1}} {_wind_fail "region must be fixed to 0 or 1"}
        dict set ctx region [dict get $region value]
    }
    set masters {}
    foreach key [list Cf $shape $view] {lappend masters [dict get $ctx $key id]}
    if {[llength [lsort -unique $masters]] != 3} {_wind_fail "Geometry coordinates must have distinct master parameters"}
    return $ctx
}

# Forward physical geometry. Angles are degrees; Cf is hemispheric coverage.
proc _wind_angles {ctx {values {}}} {
    if {$values eq ""} {
        foreach key {Cf hole skew incl view} {
            if {[dict exists $ctx $key]} {dict set values $key [lindex [_wind_out param [dict get $ctx $key id]] 0]}
        }
    }
    set cover [dict get $values Cf]
    set kind [dict get $ctx kind]
    set pi [expr {acos(-1.0)}]
    if {$kind eq "los"} {
        set s [dict get $values skew]
        set incl [dict get $values incl]
        set c [expr {$incl == 90 ? 0.0 : cos($incl*$pi/180.0)}]
        set a [expr {1.0-2.0*$s}]
        set mix [expr {(1.0-$s)*(1.0-$cover)+$s*$cover}]
        set b [expr {-($mix+$a*$c)}]
        set d [expr {sqrt(($mix-$c)**2+4.0*$s*(1.0-$s)*$c*(1.0-$c))}]
        if {$b > 0} {
            set x [expr {(-$b-$d)/(2.0*$a)}]
        } elseif {-$b+$d > 0} {
            set x [expr {2.0*(1.0-$s)*$c*(1.0-$cover)/(-$b+$d)}]
        } else {set x 0.0}
        set x [expr {max(0.0,$c-$cover,min($c,1.0-$cover,$x))}]
        set r [expr {$x+$cover}]
    } else {
        set h [dict get $values hole]
        set r [expr {1.0-$h*(1.0-$cover)}]
        set x [expr {(1.0-$h)*(1.0-$cover)}]
    }
    set tin [expr {acos(max(-1.0,min(1.0,$r)))*180.0/$pi}]
    set tout [expr {acos(max(-1.0,min(1.0,$x)))*180.0/$pi}]
    if {$kind eq "generic"} {set incl [dict get $values incl]}
    if {$kind eq "off"} {
        set v [dict get $values view]
        set incl [expr {[dict get $ctx region] == 0 ? $v*$tin : $tout+$v*(90.0-$tout)}]
    }
    return [dict create Cf $cover tin $tin tout $tout incl $incl]
}

proc windgeom {component args} {
    set ::xspec_tclout {}
    set ctx [_wind_context $component [_wind_options {*}$args]]
    set result [_wind_angles $ctx]
    dict set result uncertainty {}
    dict set result uncertainty_status unavailable
    dict set result uncertainty_state {}
    foreach quantity {Cf tin tout incl} {dict set result uncertainty_state $quantity unavailable}
    if {![catch {_wind_covariance_data} covariance]} {
        dict set result uncertainty_status local
        foreach quantity {Cf tin tout incl} {
            if {[catch {_wind_gradient $ctx $result $quantity} gradient]} {
                dict set result uncertainty_state $quantity singular
            } elseif {$gradient eq ""} {
                dict set result uncertainty_state $quantity fixed
            } elseif {![catch {_wind_covariance_error $gradient $covariance} sigma]} {
                dict set result uncertainty $quantity $sigma
                dict set result uncertainty_state $quantity local
            }
        }
    }
    puts [format "    Component %d (%s): wind geometry" $component [dict get $ctx name]]
    puts "    Quantity        Value    Unit       Covariance 1-sigma"
    foreach key {Cf tin tout incl} {
        set unit [expr {$key eq "Cf" ? "" : "deg"}]
        set error {}
        if {[dict exists $result uncertainty $key]} {set error [format %.5g [dict get $result uncertainty $key]]}
        puts [format "    %-8s %12.5f    %-3s %20s" $key [dict get $result $key] $unit $error]
    }
    puts ""
    set ::xspec_tclout $result
    return
}

# Same lower-triangle ordering/diagonal validation as the existing covariance
# procedures, extended to qualified model IDs. No dependency on xspeclab.
proc _wind_covariance_data {} {
    set ids {}
    foreach name [_wind_out modnames] {
        set prefix [expr {$name eq "_DEFAULT" ? "" : "$name:"}]
        set modelArgs [expr {$name eq "_DEFAULT" ? {} : [list $name]}]
        for {set n 1} {$n <= [_wind_out modpar {*}$modelArgs]} {incr n} {
            set id ${prefix}$n
            if {[_wind_out pfree $id] eq "T"} {lappend ids $id}
        }
    }
    set count [llength $ids]
    if {$count == 0 || $count != [_wind_out varpar]} {_wind_fail "Cannot map the free covariance parameters"}
    set lower [_wind_out covariance]
    if {[llength $lower] != $count*($count+1)/2} {_wind_fail "No current full covariance matrix"}
    foreach value $lower {
        if {![string is double -strict $value] || ![expr {abs($value)<1e300}]} {_wind_fail "Non-finite covariance"}
    }
    for {set i 0} {$i < $count} {incr i} {
        set sigma [_wind_out sigma [lindex $ids $i]]
        set diagonal [lindex $lower [expr {$i*($i+3)/2}]]
        if {$sigma < 0 || $diagonal < 0 || abs($sigma*$sigma-$diagonal)>1e-4*max(abs($diagonal),1e-30)} {
            _wind_fail "Covariance order does not match parameter sigmas"
        }
    }
    return [list $ids $lower]
}

proc _wind_covariance_error {gradient covariance} {
    lassign $covariance ids lower
    set variance 0.0
    set k 0
    for {set i 0} {$i < [llength $ids]} {incr i} {
        set left [lindex $ids $i]
        set gi [expr {[dict exists $gradient $left] ? [dict get $gradient $left] : 0.0}]
        for {set j 0} {$j <= $i} {incr j} {
            set right [lindex $ids $j]
            set gj [expr {[dict exists $gradient $right] ? [dict get $gradient $right] : 0.0}]
            set variance [expr {$variance+($i==$j ? 1.0 : 2.0)*$gi*$gj*[lindex $lower $k]}]
            incr k
        }
    }
    if {$variance < -1e-10} {_wind_fail "Negative propagated covariance variance"}
    return [expr {sqrt(max(0.0,$variance))}]
}

# Analytic angular derivatives; missing/singular gradients are not zero errors.
proc _wind_gradient {ctx geometry quantity} {
    set cover [dict get $geometry Cf]
    set kind [dict get $ctx kind]
    set shape [expr {$kind eq "los" ? "skew" : "hole"}]
    set view [expr {$kind eq "off" ? "view" : "incl"}]
    set gradient {}
    if {$quantity eq "Cf"} {set derivatives [list 1.0 0.0 0.0]} elseif {$quantity eq "incl" && $kind ne "off"} {
        set derivatives [list 0.0 0.0 1.0]
    } else {
        set pi [expr {acos(-1.0)}]
        set tin [expr {[dict get $geometry tin]*$pi/180.0}]
        set tout [expr {[dict get $geometry tout]*$pi/180.0}]
        if {abs(sin($tin))<1e-10 || abs(sin($tout))<1e-10} {_wind_fail "Singular angular covariance at a boundary"}
        set r [expr {cos($tin)}]
        set x [expr {cos($tout)}]
        set s [lindex [_wind_out param [dict get $ctx $shape id]] 0]
        if {$kind eq "los"} {
            set inclination [expr {[dict get $geometry incl]*$pi/180.0}]
            set c [expr {cos($inclination)}]
            # F=(1-s)(1-x-Cf)(c-x)-s*x*(x+Cf-c)=0.
            set fx [expr {-(1.0-$s)*(1.0+$c-$cover-2.0*$x)-$s*(2.0*$x+$cover-$c)}]
            if {abs($fx)<1e-12} {_wind_fail "Singular LoS angular derivative"}
            set dxCf [expr {((1.0-$s)*($c-$x)+$s*$x)/$fx}]
            set dxShape [expr {((1.0-$x-$cover)*($c-$x)+$x*($x+$cover-$c))/$fx}]
            set dxIncl [expr {((1.0-$s)*(1.0-$x-$cover)+$s*$x)*sin($inclination)*$pi/180.0/$fx}]
            set inner [list [expr {-180.0/$pi*($dxCf+1.0)/sin($tin)}] \
                [expr {-180.0/$pi*$dxShape/sin($tin)}] [expr {-180.0/$pi*$dxIncl/sin($tin)}]]
            set outer [list [expr {-180.0/$pi*$dxCf/sin($tout)}] \
                [expr {-180.0/$pi*$dxShape/sin($tout)}] [expr {-180.0/$pi*$dxIncl/sin($tout)}]]
        } else {
            set inner [list [expr {-180.0/$pi*$s/sin($tin)}] [expr {180.0/$pi*(1.0-$cover)/sin($tin)}] 0.0]
            set outer [list [expr {180.0/$pi*(1.0-$s)/sin($tout)}] [expr {180.0/$pi*(1.0-$cover)/sin($tout)}] 0.0]
        }
        if {$quantity eq "tin"} {set derivatives $inner} elseif {$quantity eq "tout"} {set derivatives $outer} else {
            set v [lindex [_wind_out param [dict get $ctx view id]] 0]
            set polar [expr {[dict get $ctx region]==0}]
            set base [expr {$polar ? $inner : $outer}]
            set weight [expr {$polar ? $v : 1.0-$v}]
            set derivatives [list [expr {$weight*[lindex $base 0]}] [expr {$weight*[lindex $base 1]}] \
                [expr {$polar ? [dict get $geometry tin] : 90.0-[dict get $geometry tout]}]]
        }
    }
    foreach key [list Cf $shape $view] derivative $derivatives {
        if {[dict get $ctx $key free]} {dict set gradient [dict get $ctx $key id] $derivative}
    }
    return $gradient
}

# Only a worker process changes parameters or performs fits.
proc _wind_fix {entry value} {
    set id [dict get $entry id]
    newpar $id $value
    freeze $id
}

proc _wind_constrain {ctx tin tout} {
    if {$tin > $tout} {return -code error -errorcode {WIND EXCLUDED} "tin exceeds tout"}
    set pi [expr {acos(-1.0)}]
    set r [expr {$tin == 90 ? 0.0 : cos($tin*$pi/180.0)}]
    set x [expr {$tout == 90 ? 0.0 : cos($tout*$pi/180.0)}]
    _wind_fix [dict get $ctx Cf] [expr {$r-$x}]
    if {[dict get $ctx kind] ne "los"} {
        set h [expr {$tin == 0 && $tout == 90 ? 0.5 : (1.0-$r)/(1.0-$r+$x)}]
        _wind_fix [dict get $ctx hole] $h
        return
    }
    set entry [dict get $ctx incl]
    set id [dict get $entry id]
    lassign [dict get $entry spec] value delta hardMin softMin softMax hardMax
    set lo [expr {max($tin,$hardMin)}]
    set hi [expr {min($tout,$hardMax)}]
    if {$lo > $hi || (![dict get $entry free] && ($value < $lo || $value > $hi))} {
        return -code error -errorcode {WIND EXCLUDED} "inclination constraints exclude this angular pair"
    }
    if {$tin == $tout} {
        _wind_fix [dict get $ctx skew] 0.5
        _wind_fix $entry $tin
        return
    }
    if {$lo == $hi} {
        _wind_fix $entry $lo
    } elseif {[dict get $entry free]} {
        newpar $id [expr {max($lo,min($hi,$value))}] $delta $lo $lo $hi $hi
        thaw $id
    }
    if {$tin == 0 || $tout == 90} {
        set skew [expr {$tin == 0 ? ($tout == 90 ? 0.5 : 0.0) : 1.0}]
        _wind_fix [dict get $ctx skew] $skew
    } else {
        # XSPEC integer literals denote parameters. Format all constants as reals.
        set rp [format %.17e $r]
        set xp [format %.17e $x]
        set a "(1.0-$rp)*(cosd($id)-$xp)"
        set b "$xp*($rp-cosd($id))"
        newpar [dict get $ctx skew id] = "($a)/(($a)+($b))"
    }
}

proc _wind_scan_worker {component options axes snapshot directory} {
    set result [dict create schema 1 axes $axes component $component options $options \
        statmethod [_wind_out statmethod] zerostat [_wind_out stat] \
        dof [lindex [_wind_out dof] 0] bins [lindex [_wind_out dof] 1] rows {}]
    fconfigure stdout -buffering line
    lassign $axes first second
    foreach y [dict get $second values] {
        foreach x [dict get $first values] {
            set tin [expr {[dict get $first name] eq "tin" ? $x : $y}]
            set tout [expr {[dict get $first name] eq "tout" ? $x : $y}]
            set status ok
            set message {}
            set statistic NaN
            if {$tin > $tout} {
                set status excluded
                set message "tin exceeds tout"
            } else {
                # Reload the same complete reference, including bounds and ties.
                if {[catch {
                    @$snapshot
                    query no
                    set ctx [_wind_context $component $options]
                    _wind_constrain $ctx $tin $tout
                    fit
                    set statistic [_wind_out stat]
                    if {![string is double -strict $statistic] || ![expr {abs($statistic)<1e300}]} {
                        _wind_fail "Non-finite fit statistic"
                    }
                    set actual [_wind_angles $ctx]
                    if {abs(cos([dict get $actual tin]*acos(-1.0)/180.0)-cos($tin*acos(-1.0)/180.0))>2e-9 ||
                        abs(cos([dict get $actual tout]*acos(-1.0)/180.0)-cos($tout*acos(-1.0)/180.0))>2e-9} {
                        _wind_fail "Fixed opening angles drifted during fit"
                    }
                } message details]} {
                    set statistic NaN
                    set status failed
                    if {$message eq ""} {set message "XSPEC command failed; inspect worker.log"}
                    if {[dict get $details -errorcode] eq {WIND EXCLUDED}} {set status excluded}
                } else {set message {}}
            }
            set delta NaN
            if {$status eq "ok"} {
                set delta [expr {$statistic-[dict get $result zerostat]}]
                if {$delta < -1.0 && (![dict exists $result betterfit] ||
                    $statistic < [dict get $result betterfit stat])} {
                    set candidate [file join $directory better-fit.xcm]
                    file delete -force $candidate
                    save all $candidate
                    dict set result betterfit [dict create stat $statistic tin $tin tout $tout file $candidate]
                }
            }
            dict lappend result rows [dict create values [list $x $y] stat $statistic \
                deltastat $delta status $status message $message]
            # Preserve every completed point even if a later evaluation aborts XSPEC.
            set output [open [file join $directory partial.tcl] w]
            puts $output $result
            close $output
            puts "WIND_PROGRESS [llength [dict get $result rows]]"
        }
    }
    set output [open [file join $directory result.tcl] w]
    puts $output $result
    close $output
}


# Shared full-state snapshot for scans and one-dimensional profiles.
proc _wind_snapshot {directory} {
    file mkdir $directory
    set snapshot [file join $directory snapshot.xcm]
    if {[_wind_out datasets] < 1 || [lindex [_wind_out dof] 1] < 1} {_wind_fail "Load data and responses before scanning"}
    set energyGrid [_wind_out energies 1]
    for {set spectrum 2} {$spectrum <= [_wind_out datasets]} {incr spectrum} {
        if {[_wind_out energies $spectrum] ne $energyGrid} {
            _wind_fail "Angular scans currently require a common model energy grid across spectra; set energies explicitly"
        }
    }
    set gridFile [file join $directory energies.txt]
    set output [open $gridFile w]
    puts $output [join $energyGrid \n]
    close $output
    save all $snapshot
    # save all omits an autonomous energy grid. Preserve it explicitly.
    set output [open $snapshot a]
    puts $output [list energies $gridFile]
    # Native save all formats parameter settings to six significant digits.
    # Refresh unlinked parameters with tclout's ten-digit representation.
    foreach name [_wind_out modnames] {
        set prefix [expr {$name eq "_DEFAULT" ? "" : "$name:"}]
        set modelArgs [expr {$name eq "_DEFAULT" ? {} : [list $name]}]
        for {set n 1} {$n <= [_wind_out modpar {*}$modelArgs]} {incr n} {
            set id ${prefix}$n
            if {[lindex [_wind_out plink $id] 0] eq "F"} {
                puts $output [list newpar $id {*}[_wind_out param $id]]
            }
        }
    }
    close $output
    return $snapshot
}

# Restrict one nuisance parameter while preserving its original frozen state.
proc _wind_range {entry low high} {
    set id [dict get $entry id]
    lassign [dict get $entry spec] value delta hardMin softMin softMax hardMax
    set low [expr {max($low,$hardMin)}]
    set high [expr {min($high,$hardMax)}]
    if {$low>$high || (![dict get $entry free] && ($value<$low || $value>$high))} {
        _wind_fail "Trial angle is incompatible with fixed parameters or hard limits"
    }
    if {![dict get $entry free]} {return}
    if {$low==$high} {_wind_fix $entry $low; return}
    newpar $id [expr {max($low,min($high,$value))}] $delta $low $low $high $high
    thaw $id
}

# Holding one opening angle removes Cf through the inverse geometric identity.
# All other originally free coordinates remain nuisance parameters.
proc _wind_single_constraint {ctx quantity angle} {
    set kind [dict get $ctx kind]
    if {$quantity eq "incl"} {
        if {$kind ne "off"} {_wind_fix [dict get $ctx incl] $angle; return}
        _wind_fail "Physical off-wind inclination profiling with free geometry is not supported yet"
    }
    set pi [expr {acos(-1.0)}]
    set cosine [expr {cos($angle*$pi/180.0)}]
    set cf [dict get $ctx Cf]
    if {![dict get $cf free]} {
        set cover [dict get $cf value]
        if {$quantity eq "tin"} {
            set r $cosine
            set x [expr {$r-$cover}]
        } else {set x $cosine; set r [expr {$x+$cover}]}
        if {$x<0 || $r>1} {_wind_fail "Fixed Cf excludes this trial angle"}
        if {$kind eq "los" && ![dict get $ctx skew free]} {_wind_fail "A fixed Cf and skew require a separate inclination-only profile"}
        _wind_constrain $ctx [expr {acos($r)*180.0/$pi}] [expr {acos($x)*180.0/$pi}]
        return
    }
    if {[lindex [dict get $cf spec] 2] != 0 || [lindex [dict get $cf spec] 5] != 1} {
        _wind_fail "Linked-Cf angular profiles currently require the full hard range 0..1"
    }
    set constant [format %.17e $cosine]
    # Preserve the reference's other opening angle where possible. Leaving
    # skew/hole unchanged would drive Cf toward unity near tin=0 and make a
    # physically nearby trial a numerically distant start.
    set start [_wind_angles $ctx]
    set tin [dict get $start tin]
    set tout [dict get $start tout]
    if {$quantity eq "tin"} {
        set tin $angle
        if {$tout <= $tin} {set tout 90.0}
    } else {
        set tout $angle
        if {$tin >= $tout} {set tin 0.0}
    }
    set r [expr {cos($tin*$pi/180.0)}]
    set x [expr {$tout==90 ? 0.0 : cos($tout*$pi/180.0)}]
    if {$kind eq "los"} {
        set s [dict get $ctx skew id]
        set i [dict get $ctx incl id]
        if {$quantity eq "tin"} {
            _wind_range [dict get $ctx incl] $angle 90.0
            # q=r-c; Cf=q*((1-s)*(1-r)+s*r)/((1-s)*(1-r)+s*q).
            set q "($constant-cosd($i))"
            set a "(1.0-$s)*(1.0-$constant)"
            set link "$q*(($a)+$s*$constant)/(($a)+$s*$q)"
        } else {
            _wind_range [dict get $ctx incl] 0.0 $angle
            # d=c-x; Cf=d*((1-s)*(1-x)+s*x)/((1-s)*d+s*x).
            set d "(cosd($i)-$constant)"
            set link "$d*((1.0-$s)*(1.0-$constant)+$s*$constant)/((1.0-$s)*$d+$s*$constant)"
        }
        set inclination [lindex [_wind_out param $i] 0]
        if {[dict get $ctx incl free] && ($inclination <= $tin || $inclination >= $tout)} {
            set inclination [expr {($tin+$tout)/2.0}]
            newpar $i $inclination
        }
        if {[dict get $ctx skew free]} {
            set c [expr {cos($inclination*$pi/180.0)}]
            set a [expr {(1.0-$r)*($c-$x)}]
            set b [expr {$x*($r-$c)}]
            if {$a+$b>0} {
                set spec [dict get $ctx skew spec]
                newpar $s [expr {max([lindex $spec 2],min([lindex $spec 5],$a/($a+$b)))}]
            }
        }
        if {![dict get $ctx skew free] && [dict get $ctx skew value] in {0 1 0.0 1.0}} {
            _wind_fail "A fixed endpoint skew has a singular angular chart; use the two-angle scan"
        }
    } else {
        set h [dict get $ctx hole id]
        if {$quantity eq "tin"} {
            _wind_range [dict get $ctx hole] [expr {1.0-$cosine}] 1.0
            set link "1.0-(1.0-$constant)/$h"
        } else {
            _wind_range [dict get $ctx hole] 0.0 [expr {1.0-$cosine}]
            set link "1.0-$constant/(1.0-$h)"
        }
        if {[dict get $ctx hole free]} {
            set spec [_wind_out param $h]
            set initial [expr {(1.0-$r)/(1.0-$r+$x)}]
            newpar $h [expr {max([lindex $spec 2],min([lindex $spec 5],$initial))}]
        }
    }
    newpar [dict get $cf id] = $link
}

proc _wind_profile_trial {profile angle} {
    set snapshot [dict get $profile snapshot]
    @$snapshot
    query no
    set ctx [_wind_context [dict get $profile component] [dict get $profile options]]
    set quantity [dict get $profile quantity]
    _wind_single_constraint $ctx $quantity $angle
    file delete -force $::_wine_profile_failed_trial
    save all $::_wine_profile_failed_trial
    fit
    set actual [_wind_angles $ctx]
    set pi [expr {acos(-1.0)}]
    if {abs(cos([dict get $actual $quantity]*$pi/180.0)-cos($angle*$pi/180.0))>2e-8} {
        _wind_fail "Fixed $quantity drifted at trial $angle degrees"
    }
    # A link can evaluate beyond hard limits; never interpret a clipped model
    # prediction as a valid angular constraint.
    if {[dict get $actual Cf]<-1e-9 || [dict get $actual Cf]>1.0+1e-9} {_wind_fail "Constrained Cf is outside its physical range"}
    set statistic [_wind_out stat]
    if {![string is double -strict $statistic] || ![expr {abs($statistic)<1e300}]} {_wind_fail "Non-finite profile statistic"}
    set output [open [file join [dict get $profile directory] profile-points.txt] a]
    puts $output [list $angle $statistic]
    close $output
    return $statistic
}

proc _wind_error_worker {component options quantity delta toler directory} {
    source [file join [file dirname $::_wind_source] wine_profile.tcl]
    set ::_wine_profile_output_dir $directory
    set ::_wine_profile_failed_trial [file join $directory failed-trial.xcm]
    set failed [catch {
        set original [_wind_out stat]
        fit
        set statistic [_wind_out stat]
        set ctx [_wind_context $component $options]
        set geometry [_wind_angles $ctx]
        set center [dict get $geometry $quantity]
        _wine_profile_check_improvement $center $statistic $original $toler
        set baseline [_wind_snapshot [file join $directory baseline]]
        set low 1e-5
        set high [expr {90.0-1e-5}]
        if {![dict get $ctx Cf free]} {
            set cover [dict get $ctx Cf value]
            if {$quantity eq "tin"} {set high [expr {min($high,acos($cover)*180.0/acos(-1.0))}]}
            if {$quantity eq "tout"} {set low [expr {max($low,acos(1.0-$cover)*180.0/acos(-1.0))}]}
        } elseif {[dict get $ctx kind] ne "los" && ![dict get $ctx hole free]} {
            set boundary [expr {acos(1.0-[dict get $ctx hole value])*180.0/acos(-1.0)}]
            if {$quantity eq "tin"} {set high [expr {min($high,$boundary)}]}
            if {$quantity eq "tout"} {set low [expr {max($low,$boundary)}]}
        }
        if {[dict get $ctx kind] eq "los" && ![dict get $ctx incl free]} {
            set inclination [dict get $ctx incl value]
            if {$quantity eq "tin"} {set high [expr {min($high,$inclination)}]}
            if {$quantity eq "tout"} {set low [expr {max($low,$inclination)}]}
        }
        set profile [dict create component $component options $options quantity $quantity \
            snapshot $baseline directory $directory profile_bounds [list $low $high]]
        set points [list [list $center $statistic]]
        set nativeKey {}
        if {$quantity eq "incl"} {
            set nativeKey [expr {[dict get $ctx kind] eq "off" ? "view" : "incl"}]
        } elseif {[dict get $ctx kind] eq "los" && ![dict get $ctx Cf free] && ![dict get $ctx skew free]} {
            set nativeKey incl
        }
        if {$nativeKey ne ""} {
            set id [dict get $ctx $nativeKey id]
            if {![dict get $ctx $nativeKey free]} {_wind_fail "$nativeKey is fixed"}
            error $delta $id
            _wine_profile_check_improvement $center [_wind_out stat] $statistic $toler
            lassign [_wind_out error $id] lower upper flags
            set result [dict create center $center lower $lower upper $upper flags $flags flag_format native]
            if {$quantity ne $nativeKey} {
                # Exact monotonic conversion is valid only because the other
                # geometry coordinates are fixed, not merely at their best fit.
                dict set result native_limits [list $lower $upper]
                set values {}
                foreach key {Cf hole skew incl view} {
                    if {[dict exists $ctx $key]} {dict set values $key [dict get $ctx $key value]}
                }
                foreach side {lower upper} value [list $lower $upper] {
                    if {$value != 0} {
                        dict set values $nativeKey $value
                        dict set result $side [dict get [_wind_angles $ctx $values] $quantity]
                    }
                }
                dict set result mapping monotonic-fixed-geometry
            }
        } else {
            set lower [_wine_profile_side $profile $center $statistic $delta $toler -1 points _wind_profile_trial]
            set upper [_wine_profile_side $profile $center $statistic $delta $toler 1 points _wind_profile_trial]
            set flags {}
            foreach test [list [expr {[lindex $lower 0] eq "boundary"}] [expr {[lindex $upper 0] eq "boundary"}] \
                [expr {[lindex $lower 0] eq "unstable"}] [expr {[lindex $upper 0] eq "unstable"}]] {
                append flags [expr {$test ? "T" : "F"}]
            }
            set result [dict create center $center \
                lower [expr {[lindex $lower 0] eq "boundary" ? 0 : [lindex $lower 1]}] \
                upper [expr {[lindex $upper 0] eq "boundary" ? 0 : [lindex $upper 1]}] \
                flags $flags flag_format wind sides [list $lower $upper]]
            dict set result angular_bounds [list $low $high]
        }
        dict set result status ok
        dict set result quantity $quantity
        dict set result delta $delta
        dict set result points $points
        dict set result zerostat $statistic
    } message details]
    if {$failed} {
        if {$message eq ""} {set message "XSPEC constrained fit failed; inspect worker.log"}
        set result [dict create status failed message $message code [dict get $details -errorcode]]
        if {[lrange [dict get $details -errorcode] 0 1] eq {WINE BETTER_FIT}} {
            dict set result candidate [_wind_snapshot [file join $directory better-fit]]
        }
    }
    set output [open [file join $directory result.tcl] w]
    puts $output $result
    close $output
}

proc winderror {component quantity args} {
    set ::xspec_tclout {}
    if {$quantity ni {tin tout incl}} {_wind_fail "Usage: winderror COMPONENT tin|tout|incl ?delta VALUE? ?toler VALUE? ?model NAME? ?group N?"}
    set settings [dict create delta 2.706 toler 1.0 model {} group 1]
    if {[llength $args]%2} {_wind_fail "Each winderror option requires a value"}
    foreach {key value} $args {
        if {![dict exists $settings $key]} {_wind_fail "Unknown winderror option: $key"}
        dict set settings $key $value
    }
    foreach key {delta toler} {
        set value [dict get $settings $key]
        if {![string is double -strict $value] || !($value>0 && $value<1e300)} {_wind_fail "$key must be finite and positive"}
    }
    set options [_wind_options model [dict get $settings model] group [dict get $settings group]]
    set ctx [_wind_context $component $options]
    set shape [expr {[dict get $ctx kind] eq "los" ? "skew" : "hole"}]
    if {$quantity ne "incl" && ![dict get $ctx Cf free] && [dict get $ctx Cf value]==1} {
        _wind_fail "$quantity is fixed by the frozen full-coverage geometry"
    }
    if {$quantity ne "incl" && $shape eq "hole" && ![dict get $ctx hole free] &&
        (($quantity eq "tin" && [dict get $ctx hole value]==0) ||
         ($quantity eq "tout" && [dict get $ctx hole value]==1))} {_wind_fail "$quantity is fixed by the frozen cavity coordinate"}
    if {$quantity ne "incl" && ![dict get $ctx Cf free] && ![dict get $ctx $shape free] &&
        !([dict get $ctx kind] eq "los" && [dict get $ctx incl free])} {
        _wind_fail "This profile requires at least one free coverage/shape coordinate; fixed geometry constraints are not released"
    }
    if {$quantity eq "incl" && [dict get $ctx kind] eq "off"} {
        if {[dict get $ctx Cf free] || [dict get $ctx hole free]} {
            _wind_fail "Off-wind inclination profiles require fixed Cf and hole in this version; limits on view alone are not physical inclination limits when geometry varies"
        }
        set geometry [_wind_angles $ctx]
        if {([dict get $ctx region]==0 && [dict get $geometry tin]==0) ||
            ([dict get $ctx region]==1 && [dict get $geometry tout]==90)} {_wind_fail "The selected off-wind sector has fixed physical inclination"}
    }
    set marker [file tempfile directory]
    close $marker
    file delete $directory
    file mkdir [file join $directory .xspec cache]
    set completed 0
    try {
        set snapshot [_wind_snapshot $directory]
        set script [file join $directory worker.xcm]
        set output [open $script w]
        puts $output {query no}
        puts $output [list lmod wine $::_wind_library]
        puts $output "@$snapshot"
        puts $output [list source $::_wind_source]
        puts $output [list _wind_error_worker $component $options $quantity \
            [dict get $settings delta] [dict get $settings toler] $directory]
        puts $output exit
        close $output
        puts [format "    WINE %s confidence range (delta-stat %.3g)" $quantity [dict get $settings delta]]
        set failed [catch {exec env HOME=$directory xspec - < $script > [file join $directory worker.log] 2>@1} message]
        if {$failed || ![file exists [file join $directory result.tcl]]} {_wind_fail "Profile stopped. Diagnostics: $directory"}
        set input [open [file join $directory result.tcl] r]
        set result [string trim [read $input]]
        close $input
        if {[dict get $result status] ne "ok"} {
            if {[lrange [dict get $result code] 0 1] eq {WINE BETTER_FIT}} {
                _wind_fail "A better minimum exceeds toler. Inspect [dict get $result candidate] and refit the starting solution before reporting errors"
            }
            _wind_fail "[dict get $result message]. Diagnostics: $directory"
        }
        puts [format "    Reference %s (component %d): %.2f deg" $quantity $component [dict get $result center]]
        puts "    Component  Quantity       Lower        Upper       Flags"
        puts [format "    %9d  %-8s %12.2f %12.2f       %s" $component $quantity \
            [dict get $result lower] [dict get $result upper] [dict get $result flags]]
        if {[dict get $result flag_format] eq "wind" && [string match *T* [dict get $result flags]]} {
            puts "    Warning: check sides for pegged or numerically unstable limits."
            puts "    A boundary value displayed as 0 is not a measured confidence bound."
        } elseif {[dict get $result flag_format] eq "native" && [string match *T* [dict get $result flags]]} {
            puts "    Warning: native XSPEC error flags are set; inspect flags/native_limits."
        }
        set completed 1
        puts ""
        set ::xspec_tclout $result
        return
    } finally {
        if {$completed} {file delete -force $directory}
    }
}

proc windsteppar {component args} {
    set ::xspec_tclout {}
    if {[llength $args] < 8} {_wind_fail "Usage: windsteppar COMPONENT tin MIN MAX INTERVALS tout MIN MAX INTERVALS ?model NAME? ?group N?"}
    set axes {}
    foreach offset {0 4} {
        lassign [lrange $args $offset [expr {$offset+3}]] name lo hi steps
        if {$name ni {tin tout} || ![string is double -strict $lo] || ![string is double -strict $hi] ||
            ![string is integer -strict $steps] || $steps<1 || !($lo>=0 && $hi<=90 && $lo<$hi)} {
            _wind_fail "Each axis requires tin/tout, 0 <= MIN < MAX <= 90, and positive integer INTERVALS"
        }
        set values {}
        for {set i 0} {$i <= $steps} {incr i} {lappend values [expr {$lo+($hi-$lo)*$i/double($steps)}]}
        lappend axes [dict create name $name unit deg values $values]
    }
    if {[dict get [lindex $axes 0] name] eq [dict get [lindex $axes 1] name]} {_wind_fail "Specify both tin and tout"}
    set options [_wind_options {*}[lrange $args 8 end]]
    _wind_context $component $options
    set marker [file tempfile directory]
    close $marker
    file delete $directory
    file mkdir [file join $directory .xspec cache]
    set snapshot [file join $directory snapshot.xcm]
    set script [file join $directory worker.xcm]
    set success 0
    try {
        _wind_snapshot $directory
        set output [open $script w]
        puts $output {query no}
        puts $output [list lmod wine $::_wind_library]
        puts $output "@$snapshot"
        puts $output [list source $::_wind_source]
        puts $output [list _wind_scan_worker $component $options $axes $snapshot $directory]
        puts $output exit
        close $output
        puts "    WINE angular scan: fitting in an isolated XSPEC session"
        set total [expr {[llength [dict get [lindex $axes 0] values]] *
            [llength [dict get [lindex $axes 1] values]]}]
        set log [open [file join $directory worker.log] w]
        set failed [catch {
            set worker [open [list | env HOME=$directory xspec - < $script 2>@1] r]
            try {
                while {[gets $worker line] >= 0} {
                    puts $log $line
                    if {[regexp {^WIND_PROGRESS ([0-9]+)$} $line -> done]} {
                        puts [format "    Completed %d/%d points (%3.0f%%)" $done $total \
                            [expr {100.0*$done/$total}]]
                    }
                }
            } finally {close $worker}
        } message]
        close $log
        puts ""
        if {$failed || ![file exists [file join $directory result.tcl]]} {
            _wind_fail "Angular scan stopped ($message). Diagnostics and completed points: $directory"
        }
        set input [open [file join $directory result.tcl] r]
        set result [string trim [read $input]]
        close $input
        puts "    tin       tout          Statistic      Status"
        set failedPoints 0
        foreach row [dict get $result rows] {
            lassign [dict get $row values] x y
            set tin [expr {[dict get [lindex $axes 0] name] eq "tin" ? $x : $y}]
            set tout [expr {[dict get [lindex $axes 0] name] eq "tout" ? $x : $y}]
            puts [format "    %7.2f %7.2f %18s      %s" $tin $tout \
                [expr {[dict get $row status] eq "ok" ? [format %.8g [dict get $row stat]] : "---"}] \
                [dict get $row status]]
            if {[dict get $row status] eq "failed"} {incr failedPoints}
        }
        if {$failedPoints || [dict exists $result betterfit]} {
            dict set result diagnostics $directory
            if {$failedPoints} {puts "    Warning: $failedPoints fits failed; inspect $directory"}
            if {[dict exists $result betterfit]} {
                puts "    Warning: a fit improves the reference statistic by more than 1.0."
                puts "    Inspect [dict get $result betterfit file] before interpreting confidence contours."
            }
        } else {set success 1}
        puts ""
        set ::xspec_tclout $result
        return
    } finally {
        if {$success} {file delete -force $directory}
    }
}
