# WINE viewing-angle covariance error, using XSPEC's fit covariance directly.

# Replace a geometric parameter with its independent XSPEC parameter, or its
# fixed numeric value. Follow simple identity ties between data groups.
proc _wine_cov_operand {number {visited {}}} {
    if {[lsearch -exact $visited $number] >= 0} {
        return -code error "Cyclic parameter link at $number"
    }
    tclout plink $number
    set link $xspec_tclout
    if {[lindex $link 0] eq "T"} {
        if {![regexp {^T\s+=\s+p([0-9]+)$} $link -> target]} {
            return -code error "Parameter $number has a non-identity geometric link ($link)"
        }
        return [_wine_cov_operand $target [linsert $visited end $number]]
    }
    tclout param $number
    set value [lindex $xspec_tclout 0]
    if {![string is double -strict $value]} {
        return -code error "Parameter $number has no numeric value"
    }
    tclout pfree $number
    set id [expr {$xspec_tclout eq "T" ? $number : ""}]
    return [list $value $id]
}

# Keep the same XSPEC free-parameter order as the lower covariance triangle.
proc _wine_covariance_data {} {
    tclout varpar
    set nfree $xspec_tclout
    tclout modpar
    set npar $xspec_tclout
    set ids {}
    for {set par 1} {$par <= $npar} {incr par} {
        tclout pfree $par
        if {$xspec_tclout ne "T"} {continue}
        tclout plink $par
        if {[lindex $xspec_tclout 0] eq "F"} {lappend ids $par}
    }
    if {[llength $ids] != $nfree} {
        return -code error "Cannot map all XSPEC variable parameters to model parameter numbers"
    }
    if {$nfree == 0 || [catch {tclout covariance}]} {
        return -code error "No covariance is available from the most recent fit"
    }
    set lower $xspec_tclout
    if {[llength $lower] != $nfree * ($nfree + 1) / 2} {
        return -code error "XSPEC covariance has [llength $lower] entries; expected [expr {$nfree * ($nfree + 1) / 2}]"
    }
    foreach value $lower {
        if {![string is double -strict $value] || ![expr {abs($value) < 1.0e300}]} {
            return -code error "XSPEC covariance contains a non-finite value: $value"
        }
    }
    for {set i 0} {$i < $nfree} {incr i} {
        tclout sigma [lindex $ids $i]
        set sigma $xspec_tclout
        set diagonal [lindex $lower [expr {$i * ($i + 3) / 2}]]
        if {$sigma < 0 || abs($sigma * $sigma - $diagonal) >
            1.0e-4 * max(abs($diagonal), 1.0e-30)} {
            return -code error "XSPEC covariance order does not match parameter [lindex $ids $i]"
        }
    }
    return [list $ids $lower]
}

proc _wine_cov_gradient {gradientVar operand derivative} {
    upvar 1 $gradientVar gradient
    set id [lindex $operand 1]
    if {$id eq ""} {return}
    if {![dict exists $gradient $id]} {dict set gradient $id 0.0}
    dict set gradient $id [expr {[dict get $gradient $id] + $derivative}]
}

# Return the local 1-sigma covariance error on physical inclination in degrees.
# Usage: wineinclcov N, with N an incl, uincl, or uview parameter number.
# Print the angle and error; keep full error precision in $xspec_tclout.
proc wineinclcov {number} {
    if {![string is integer -strict $number] || $number < 1} {
        return -code error "Usage: wineinclcov <incl|uincl|uview XSPEC parameter number>"
    }
    tclout pinfo $number
    set name [lindex $xspec_tclout 0]
    if {$name ni {incl uincl uview}} {
        return -code error "Parameter $number is $name, expected a WINE incl, uincl, or uview"
    }
    set view [_wine_cov_operand $number]
    set u [lindex $view 0]
    set gradient {}
    if {$name eq "incl"} {
        set value $u
        _wine_cov_gradient gradient $view 1.0
    } else {
        tclout pinfo [expr {$number + 1}]
        if {[lindex $xspec_tclout 0] ne "tout"} {
            return -code error "Expected tout immediately after $name"
        }
        tclout pinfo [expr {$number + 2}]
        if {[lindex $xspec_tclout 0] ne "tin"} {
            return -code error "Expected tin two parameters after $name"
        }
        set tout [_wine_cov_operand [expr {$number + 1}]]
        set tin [_wine_cov_operand [expr {$number + 2}]]
        set outer [lindex $tout 0]
        set inner [lindex $tin 0]
        if {$name eq "uincl"} {
            set value [expr {$inner + $u * ($outer - $inner)}]
            _wine_cov_gradient gradient $view [expr {$outer - $inner}]
            _wine_cov_gradient gradient $tout $u
            _wine_cov_gradient gradient $tin [expr {1.0 - $u}]
        } else {
            tclout pinfo [expr {$number + 3}]
            if {[lindex $xspec_tclout 0] ne "region"} {
                return -code error "Expected region three parameters after uview"
            }
            tclout param [expr {$number + 3}]
            set region [lindex $xspec_tclout 0]
            if {$region == 0} {
                set value [expr {(1.0 - $u) * $inner}]
                _wine_cov_gradient gradient $view [expr {-$inner}]
                _wine_cov_gradient gradient $tin [expr {1.0 - $u}]
            } elseif {$region == 1} {
                set value [expr {$outer + $u * (90.0 - $outer)}]
                _wine_cov_gradient gradient $view [expr {90.0 - $outer}]
                _wine_cov_gradient gradient $tout [expr {1.0 - $u}]
            } else {
                return -code error "WINE region must be 0 or 1"
            }
        }
    }
    lassign [_wine_covariance_data] ids lower
    set variance 0.0
    set scale 0.0
    dict for {par_i derivative_i} $gradient {
        set i [lsearch -exact $ids $par_i]
        if {$i < 0} {return -code error "Parameter $par_i is absent from the fit covariance"}
        dict for {par_j derivative_j} $gradient {
            set j [lsearch -exact $ids $par_j]
            if {$j < 0} {return -code error "Parameter $par_j is absent from the fit covariance"}
            set row [expr {max($i, $j)}]
            set col [expr {min($i, $j)}]
            set cov [lindex $lower [expr {$row * ($row + 1) / 2 + $col}]]
            set term [expr {$derivative_i * $cov * $derivative_j}]
            set variance [expr {$variance + $term}]
            set scale [expr {$scale + abs($term)}]
        }
    }
    if {$variance < -1.0e-10 * $scale} {
        return -code error "Propagated variance is negative: $variance"
    }
    set sigma [expr {sqrt(max(0.0, $variance))}]
    set ::xspec_tclout $sigma
    puts [format "    Covariance: incl = %.2f +/- %.2f deg" $value $sigma]
    return $::xspec_tclout
}

# Convert native XSPEC error limits to physical viewing-angle degrees.
# Usage: wineinclerror N ?delta-stat? (default XSPEC 90% one-parameter: 2.706).
# This is a physical-angle profile only for direct incl or frozen, unlinked
# opening angles. $xspec_tclout retains full-precision limits and status.
proc wineinclerror {number {delta 2.706}} {
    if {![string is integer -strict $number] || $number < 1 ||
        ![string is double -strict $delta] || $delta <= 0} {
        return -code error "Usage: wineinclerror <incl|uincl|uview XSPEC parameter number> ?delta-stat?"
    }
    tclout pinfo $number
    set name [lindex $xspec_tclout 0]
    if {$name ni {incl uincl uview}} {
        return -code error "Parameter $number is $name, expected incl, uincl, or uview"
    }
    tclout pfree $number
    if {$xspec_tclout ne "T"} {
        return -code error "Parameter $number must be free and independent for XSPEC error"
    }
    set outer 0.0
    set inner 0.0
    set region 0
    if {$name ne "incl"} {
        foreach {offset expected} {1 tout 2 tin} {
            set id [expr {$number + $offset}]
            tclout pinfo $id
            if {[lindex $xspec_tclout 0] ne $expected} {
                return -code error "Expected $expected at XSPEC parameter $id"
            }
            tclout plink $id
            if {[lindex $xspec_tclout 0] eq "T"} {
                return -code error "Opening angles must be fixed for wineinclerror; parameter $id is linked"
            }
            tclout pfree $id
            if {$xspec_tclout ne "F"} {
                return -code error "Opening angles must be frozen for wineinclerror; parameter $id is free"
            }
            tclout param $id
            if {$expected eq "tout"} {set outer [lindex $xspec_tclout 0]}
            if {$expected eq "tin"} {set inner [lindex $xspec_tclout 0]}
        }
        if {$name eq "uview"} {
            set id [expr {$number + 3}]
            tclout pinfo $id
            if {[lindex $xspec_tclout 0] ne "region"} {
                return -code error "Expected region at XSPEC parameter $id"
            }
            tclout param $id
            set region [lindex $xspec_tclout 0]
            if {$region != 0 && $region != 1} {
                return -code error "WINE region must be 0 or 1"
            }
        }
        if {$outer <= $inner} {
            return -code error "WINE geometry requires tout > tin"
        }
    }
    # In XSPEC, `error` is the native fitting command, not Tcl's error routine.
    if {[catch {error $delta $number}]} {
        return -code error "XSPEC error failed for parameter $number; inspect the XSPEC messages above"
    }
    if {[catch {tclout error $number}]} {
        return -code error "XSPEC returned no profile limits for parameter $number"
    }
    lassign $xspec_tclout low high status
    if {![string is double -strict $low] || ![string is double -strict $high]} {
        return -code error "XSPEC did not return numeric error bounds for parameter $number: $xspec_tclout"
    }
    if {$name eq "uincl"} {
        set limits [list [expr {$inner + $low * ($outer - $inner)}] \
                         [expr {$inner + $high * ($outer - $inner)}]]
    } elseif {$name eq "uview" && $region == 0} {
        set limits [list [expr {(1.0 - $high) * $inner}] \
                         [expr {(1.0 - $low) * $inner}]]
    } elseif {$name eq "uview"} {
        set limits [list [expr {$outer + $low * (90.0 - $outer)}] \
                         [expr {$outer + $high * (90.0 - $outer)}]]
    } else {
        set limits [list $low $high]
    }
    set ::xspec_tclout [concat $limits [list $status]]
    puts [format "    Viewing-angle limits (deg; delta-stat %.3g): %.2f, %.2f (%s)" \
        $delta [lindex $limits 0] [lindex $limits 1] $status]
    return $::xspec_tclout
}
