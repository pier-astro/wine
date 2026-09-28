# Profile a WINE viewing angle in an isolated XSPEC process. The calling
# session is never fitted, relinked, or assigned new parameter limits.

set _wine_profile_path [file normalize [info script]]
if {[file type $_wine_profile_path] eq "link"} {
    set _wine_profile_path [file normalize [file join \
        [file dirname $_wine_profile_path] [file readlink $_wine_profile_path]]]
}
set ::_wine_profile_source $_wine_profile_path
set ::_wine_profile_xspec [file normalize [file join \
    [file dirname $_wine_profile_path] .. xspec]]
unset _wine_profile_path

proc _wine_profile_fail {message} {return -code error $message}

proc _wine_profile_master {number {seen {}}} {
    if {[lsearch -exact $seen $number] >= 0} {
        _wine_profile_fail "cyclic link involving parameter $number"
    }
    tclout plink $number
    if {[lindex $xspec_tclout 0] ne "T"} {return $number}
    if {![regexp {^T\s+=\s+p([0-9]+)$} $xspec_tclout -> master]} {
        _wine_profile_fail "parameter $number has a non-identity link: $xspec_tclout"
    }
    return [_wine_profile_master $master [linsert $seen end $number]]
}

proc _wine_profile_parameter {number expected} {
    tclout pinfo $number
    if {[lindex $xspec_tclout 0] ne $expected} {
        _wine_profile_fail "expected $expected at parameter $number"
    }
    set master [_wine_profile_master $number]
    tclout pinfo $master
    if {[lindex $xspec_tclout 0] ne $expected} {
        _wine_profile_fail "parameter $number is linked to a different kind of parameter ($master)"
    }
    tclout param $master
    set spec $xspec_tclout
    if {[llength $spec] != 6} {
        _wine_profile_fail "cannot read six XSPEC settings for parameter $master"
    }
    tclout pfree $master
    return [dict create id $master spec $spec free [expr {$xspec_tclout eq "T"}]]
}

proc _wine_profile_context {number} {
    if {![string is integer -strict $number] || $number < 1} {
        _wine_profile_fail "usage: wineinclprofile ?delta-stat? <uincl|uview parameter number>"
    }
    tclout pinfo $number
    set kind [lindex $xspec_tclout 0]
    if {$kind ni {uincl uview}} {
        _wine_profile_fail "parameter $number must be a WINE uincl or uview"
    }
    set u [_wine_profile_parameter $number $kind]
    if {![dict get $u free]} {
        _wine_profile_fail "$kind must be free (or identity-linked to a free $kind)"
    }
    set tout [_wine_profile_parameter [expr {$number + 1}] tout]
    set tin [_wine_profile_parameter [expr {$number + 2}] tin]
    set uspec [dict get $u spec]
    if {[lindex $uspec 2] > 0.0 || [lindex $uspec 5] < 1.0} {
        _wine_profile_fail "$kind must retain its full physical hard range 0..1"
    }
    set region -1
    if {$kind eq "uview"} {
        tclout pinfo [expr {$number + 3}]
        if {[lindex $xspec_tclout 0] ne "region"} {
            _wine_profile_fail "expected region three parameters after uview"
        }
        tclout param [expr {$number + 3}]
        set region [lindex $xspec_tclout 0]
        if {$region ni {0 1}} {_wine_profile_fail "WINE region must be 0 or 1"}
        tclout pfree [expr {$number + 3}]
        if {$xspec_tclout ne "F"} {_wine_profile_fail "WINE region must be frozen"}
    }
    set ids [list [dict get $u id] [dict get $tout id] [dict get $tin id]]
    if {[llength [lsort -unique $ids]] != 3} {
        _wine_profile_fail "view and opening angles must have different free masters"
    }
    return [dict create kind $kind region $region u $u tout $tout tin $tin]
}

proc _wine_profile_angle {ctx} {
    foreach name {u tout tin} {
        set id [dict get $ctx $name id]
        tclout param $id
        set $name [lindex $xspec_tclout 0]
    }
    set kind [dict get $ctx kind]
    if {$kind eq "uincl"} {return [expr {$tin + $u * ($tout - $tin)}]}
    if {[dict get $ctx region] == 0} {return [expr {(1.0 - $u) * $tin}]}
    return [expr {$tout + $u * (90.0 - $tout)}]
}

# Set an angle to value with the baseline bounds, then optionally restrict one
# hard limit to the trial inclination. Only the worker process calls this.
proc _wine_profile_set_angle {entry value {newMin {}} {newMax {}}} {
    set id [dict get $entry id]
    lassign [dict get $entry spec] ignored delta hardMin softMin softMax hardMax
    if {$newMin ne ""} {set hardMin [expr {max($hardMin, $newMin)}]}
    if {$newMax ne ""} {set hardMax [expr {min($hardMax, $newMax)}]}
    if {$hardMin >= $hardMax || $value < $hardMin || $value > $hardMax} {
        _wine_profile_fail "no allowed range for geometric parameter $id at this inclination"
    }
    set softMin [expr {max($hardMin, min($softMin, $value))}]
    set softMax [expr {min($hardMax, max($softMax, $value))}]
    newpar $id $value $delta $hardMin $softMin $softMax $hardMax
}

proc _wine_profile_trial {ctx inclination} {
    set uentry [dict get $ctx u]
    set tinEntry [dict get $ctx tin]
    set toutEntry [dict get $ctx tout]
    set uid [dict get $uentry id]
    set tinid [dict get $tinEntry id]
    set toutid [dict get $toutEntry id]

    # Remove the previous trial's link and restore the original hard bounds.
    tclout plink $uid
    if {[lindex $xspec_tclout 0] eq "T"} {untie $uid}
    newpar $uid {*}[dict get $uentry spec]
    foreach name {tin tout} {
        set entry [dict get $ctx $name]
        tclout param [dict get $entry id]
        _wine_profile_set_angle $entry [lindex $xspec_tclout 0]
    }
    tclout param $tinid
    set tin [lindex $xspec_tclout 0]
    tclout param $toutid
    set tout [lindex $xspec_tclout 0]
    set kind [dict get $ctx kind]
    if {$kind eq "uincl"} {
        set targetTin [expr {min($tin, $inclination)}]
        set targetTout [expr {max($tout, $inclination)}]
        if {(![dict get $tinEntry free] && $targetTin != $tin) ||
            (![dict get $toutEntry free] && $targetTout != $tout)} {
            _wine_profile_fail "trial inclination lies outside frozen cone angles"
        }
        if {$targetTout > $tout} {_wine_profile_set_angle $toutEntry $targetTout}
        if {$targetTin < $tin} {_wine_profile_set_angle $tinEntry $targetTin}
        _wine_profile_set_angle $tinEntry $targetTin {} $inclination
        _wine_profile_set_angle $toutEntry $targetTout $inclination
        set link [format "(%.16g-p%d)/(p%d-p%d)" $inclination $tinid $toutid $tinid]
    } elseif {[dict get $ctx region] == 0} {
        set targetTin [expr {max($tin, $inclination)}]
        set targetTout [expr {max($tout, $targetTin + 1.0e-5)}]
        if {(![dict get $tinEntry free] && $targetTin != $tin) ||
            (![dict get $toutEntry free] && $targetTout != $tout)} {
            _wine_profile_fail "trial inclination lies outside frozen polar-cavity angles"
        }
        if {$targetTout > $tout} {_wine_profile_set_angle $toutEntry $targetTout}
        if {$targetTin > $tin} {_wine_profile_set_angle $tinEntry $targetTin}
        _wine_profile_set_angle $tinEntry $targetTin $inclination
        _wine_profile_set_angle $toutEntry $targetTout
        set link [format "1.0-%.16g/p%d" $inclination $tinid]
    } else {
        set targetTout [expr {min($tout, $inclination)}]
        set targetTin [expr {min($tin, $targetTout - 1.0e-5)}]
        if {(![dict get $tinEntry free] && $targetTin != $tin) ||
            (![dict get $toutEntry free] && $targetTout != $tout)} {
            _wine_profile_fail "trial inclination lies outside frozen equatorial angles"
        }
        if {$targetTin < $tin} {_wine_profile_set_angle $tinEntry $targetTin}
        if {$targetTout < $tout} {_wine_profile_set_angle $toutEntry $targetTout}
        _wine_profile_set_angle $toutEntry $targetTout {} $inclination
        _wine_profile_set_angle $tinEntry $targetTin
        set link [format "(%.16g-p%d)/(90.0-p%d)" $inclination $toutid $toutid]
    }
    newpar $uid = $link
    # Save the exact constrained state before XSPEC can abort inside a fit.
    file delete -force $::_wine_profile_failed_trial
    save all $::_wine_profile_failed_trial
    if {[catch {fit} message]} {
        _wine_profile_fail "XSPEC trial fit failed at inclination $inclination: $message"
    }
    tclout stat
    set statistic $xspec_tclout
    set actual [_wine_profile_angle $ctx]
    if {abs($actual - $inclination) > 1.0e-5} {
        _wine_profile_fail "fixed viewing angle drifted to $actual instead of $inclination"
    }
    tclout param $uid
    set u [lindex $xspec_tclout 0]
    tclout param $tinid
    set finalTin [lindex $xspec_tclout 0]
    tclout param $toutid
    set finalTout [lindex $xspec_tclout 0]
    if {$u < -1.0e-7 || $u > 1.0 + 1.0e-7 ||
        $finalTout <= $finalTin ||
        ($kind eq "uincl" &&
         ($inclination < $finalTin - 1.0e-6 ||
          $inclination > $finalTout + 1.0e-6)) ||
        ($kind eq "uview" && [dict get $ctx region] == 0 &&
         $inclination > $finalTin + 1.0e-6) ||
        ($kind eq "uview" && [dict get $ctx region] == 1 &&
         $inclination < $finalTout - 1.0e-6)} {
        _wine_profile_fail "trial fit ended outside its physical viewing geometry at inclination $inclination"
    }
    return $statistic
}

proc _wine_profile_check_improvement {angle stat bestStat toler} {
    set gain [expr {$bestStat - $stat}]
    if {$gain <= $toler} {return}
    set index 1
    while {1} {
        set name [expr {$index == 1 ? "better-fit.xcm" : "better-fit-$index.xcm"}]
        set candidate [file join $::_wine_profile_output_dir $name]
        if {![file exists $candidate]} {break}
        incr index
    }
    if {[catch {save all $candidate} message]} {
        _wine_profile_fail "better fit found, but could not save $candidate: $message"
    }
    return -code error -errorcode [list WINE BETTER_FIT $angle $stat $gain $toler $name] \
        "a better fit was found"
}

proc _wine_profile_side {ctx center bestStat delta toler direction pointsVar} {
    upvar 1 $pointsVar points
    set previous $center
    set previousStat $bestStat
    set step 1.0
    set low 1.0e-5
    set high [expr {90.0 - 1.0e-5}]
    foreach name {tin tout} {
        set entry [dict get $ctx $name]
        if {[dict get $entry free]} {continue}
        set value [lindex [dict get $entry spec] 0]
        set kind [dict get $ctx kind]
        if {$kind eq "uincl"} {
            if {$name eq "tin"} {set low [expr {max($low, $value)}]}
            if {$name eq "tout"} {set high [expr {min($high, $value)}]}
        } elseif {[dict get $ctx region] == 0 && $name eq "tin"} {
            set high [expr {min($high, $value)}]
        } elseif {[dict get $ctx region] == 1 && $name eq "tout"} {
            set low [expr {max($low, $value)}]
        }
    }
    set edge [expr {$direction < 0 ? $low : $high}]
    for {set trial 0} {$trial < 14} {incr trial} {
        set angle [expr {$direction < 0 ? max($edge, $center - $step) : min($edge, $center + $step)}]
        if {abs($angle - $previous) < 1.0e-8} {return [list boundary $previous]}
        set stat [_wine_profile_trial $ctx $angle]
        lappend points [list $angle $stat]
        set excess [expr {$stat - $bestStat}]
        _wine_profile_check_improvement $angle $stat $bestStat $toler
        if {$excess >= $delta} {
            set inside $previous
            set outside $angle
            set insideStat $previousStat
            set outsideStat $stat
            for {set refine 0} {$refine < 22} {incr refine} {
                set middle [expr {($inside + $outside) / 2.0}]
                set stat [_wine_profile_trial $ctx $middle]
                lappend points [list $middle $stat]
                set excess [expr {$stat - $bestStat}]
                _wine_profile_check_improvement $middle $stat $bestStat $toler
                if {$excess >= $delta} {
                    set outside $middle
                    set outsideStat $stat
                } else {
                    set inside $middle
                    set insideStat $stat
                }
                if {abs($outside - $inside) < 0.01} {
                    if {abs($excess - $delta) <= 0.01} {
                        return [list ok $middle]
                    }
                    set mismatch [expr {max(abs($insideStat - $bestStat - $delta), \
                        abs($outsideStat - $bestStat - $delta))}]
                    if {$mismatch <= $toler} {
                        return [list approximate [expr {($inside + $outside) / 2.0}] \
                            [list $inside $outside] [list $insideStat $outsideStat]]
                    }
                }
            }
            if {abs($outside - $inside) < 0.001} {
                set mismatch [expr {max(abs($insideStat - $bestStat - $delta), \
                    abs($outsideStat - $bestStat - $delta))}]
                set status [expr {$mismatch <= $toler ? "approximate" : "unstable"}]
                return [list $status [expr {($inside + $outside) / 2.0}] \
                    [list $inside $outside] [list $insideStat $outsideStat]]
            }
            _wine_profile_fail "profile search did not converge on the [expr {$direction < 0 ? "lower" : "upper"}] side"
        }
        if {$angle == $edge} {return [list boundary $edge]}
        set previous $angle
        set previousStat $stat
        set step [expr {$step * 2.0}]
    }
    _wine_profile_fail "could not bracket the profile limit"
}

proc _wine_profile_compute {number delta toler} {
    set ctx [_wine_profile_context $number]
    if {[catch {fit} message]} {_wine_profile_fail "XSPEC baseline fit failed: $message"}
    tclout stat
    set bestStat $xspec_tclout
    set center [_wine_profile_angle $ctx]
    if {$center <= 0.0 || $center >= 90.0} {
        _wine_profile_fail "best inclination is on the 0/90-degree boundary"
    }
    set points [list [list $center $bestStat]]
    set lower [_wine_profile_side $ctx $center $bestStat $delta $toler -1 points]
    set upper [_wine_profile_side $ctx $center $bestStat $delta $toler 1 points]
    return [list $center $bestStat $lower $upper $points]
}

proc _wine_profile_run {number delta resultFile outputDir {trialAngle {}} {fitDelta 0.001} {toler 1.0}} {
    set ::_wine_profile_failed_trial [file join [file dirname $resultFile] failed-trial.xcm]
    set ::_wine_profile_output_dir $outputDir
    xset delta $fitDelta
    if {$trialAngle eq ""} {
        set failed [catch {_wine_profile_compute $number $delta $toler} result options]
    } else {
        set failed [catch {
            set ctx [_wine_profile_context $number]
            _wine_profile_trial $ctx $trialAngle
        } result options]
    }
    if {$failed} {
        set payload [list failed $result [dict get $options -errorcode]]
        puts stderr "WINE profile diagnostic: [dict get $options -errorinfo]"
    } elseif {$trialAngle ne ""} {
        set payload [list ok $result]
    } else {
        set payload [linsert $result 0 ok]
    }
    set output [open $resultFile w]
    puts $output $payload
    close $output
}

# Parse the public XSPEC-like syntax. A leading real is delta fit statistic;
# `wineinclprofile 1. 2` requests delta-stat 1 for parameter 2.
# `toler` is the absolute statistic tolerance for quiet approximate limits
# and for reporting a newly found minimum; it does not change delta-stat.
proc _wine_profile_arguments {args} {
    set fitDelta 0.001
    set toler 1.0
    set seen [dict create]
    while {[llength $args] >= 2 && [lindex $args end-1] in {-fit-delta toler}} {
        set option [lindex $args end-1]
        set value [lindex $args end]
        set args [lrange $args 0 end-2]
        if {[dict exists $seen $option]} {
            _wine_profile_fail "duplicate $option option"
        }
        dict set seen $option 1
        if {![string is double -strict $value] || $value != $value ||
            $value <= 0 || $value >= 1e300} {
            _wine_profile_fail "$option must be a positive finite number"
        }
        if {$option eq "-fit-delta"} {
            set fitDelta $value
        } else {
            set toler $value
        }
    }
    if {$fitDelta >= 1} {
        _wine_profile_fail "-fit-delta must be between 0 and 1"
    }
    set trialAngle {}
    if {[llength $args] == 3 && [lindex $args 0] eq "-trial"} {
        lassign [lrange $args 1 end] trialAngle number
        if {![string is double -strict $trialAngle] ||
            $trialAngle != $trialAngle ||
            $trialAngle <= 0 || $trialAngle >= 90 ||
            ![string is integer -strict $number] || $number < 1} {
            _wine_profile_fail "usage: wineinclprofile -trial <angle-deg> <uincl|uview parameter number>"
        }
        return [list $number {} $trialAngle $fitDelta $toler]
    } elseif {[llength $args] == 1} {
        set number [lindex $args 0]
        set delta 2.706
    } elseif {[llength $args] == 2} {
        lassign $args delta number
        if {[string is integer -strict $delta]} {
            _wine_profile_fail "write delta-stat as a real number, for example 1.0"
        }
    } else {
        _wine_profile_fail "usage: wineinclprofile ?delta-stat? <uincl|uview parameter number>"
    }
    if {![string is integer -strict $number] || $number < 1 ||
        ![string is double -strict $delta] || $delta != $delta ||
        $delta <= 0 || $delta >= 1e300} {
        _wine_profile_fail "usage: wineinclprofile ?delta-stat? <uincl|uview parameter number>"
    }
    return [list $number $delta $trialAngle $fitDelta $toler]
}

# Display XSPEC-style limits and expose full-precision values for Tcl scripts.
# Flags are lower boundary, upper boundary, lower instability, upper instability.
# An approximate crossing within toler carries no warning flag.
proc _wine_profile_present {number delta center lower upper} {
    set lowerStatus [lindex $lower 0]
    set upperStatus [lindex $upper 0]
    set lowerValue [expr {$lowerStatus eq "boundary" ? 0 : [lindex $lower 1]}]
    set upperValue [expr {$upperStatus eq "boundary" ? 0 : [lindex $upper 1]}]
    set flags {}
    foreach marked [list \
        [expr {$lowerStatus eq "boundary"}] [expr {$upperStatus eq "boundary"}] \
        [expr {$lowerStatus eq "unstable"}] [expr {$upperStatus eq "unstable"}]] {
        append flags [expr {$marked ? "T" : "F"}]
    }
    set ::xspec_tclout [list $lowerValue $upperValue $flags]
    puts [format " Reference Inclination (p%d): %.2f deg" $number $center]
    puts [format " Parameter     Confidence Range (%.3g)" $delta]
    set shownLower [expr {$lowerStatus eq "boundary" ? "0" : [format %.2f $lowerValue]}]
    set shownUpper [expr {$upperStatus eq "boundary" ? "0" : [format %.2f $upperValue]}]
    puts [format "    %5d %12s %12s    (%.2f,%.2f)" \
        $number $shownLower $shownUpper \
        [expr {$lowerValue - $center}] [expr {$upperValue - $center}]]
    foreach side {lower upper} status [list $lowerStatus $upperStatus] \
            limit [list $lower $upper] {
        if {$status eq "boundary"} {
            puts [format "***Warning: p%d reached the %s allowed limit; no confidence bound on this side." \
                $number $side]
        } elseif {$status eq "unstable"} {
            puts [format "***Warning: %s limit is numerically unstable; check wineinclprofile -trial %.2f %d." \
                $side [lindex $limit 1] $number]
        }
    }
}

# Profile the physical WINE viewing angle using constrained XSPEC fits.
# Usage: wineinclprofile ?delta-stat? N ?toler value? ?-fit-delta fraction?
# Default delta-stat is 2.706; N is a uincl or uview parameter number.
# Use `-trial angle N` for one diagnostic fit. The calling session is preserved.
# $xspec_tclout holds full-precision limits and a four-letter diagnostic flag.
proc wineinclprofile {args} {
    lassign [_wine_profile_arguments {*}$args] number delta trialAngle fitDelta toler
    set marker [file tempfile temporary]
    close $marker
    file delete $temporary
    file mkdir $temporary
    set snapshot [file join $temporary snapshot.xcm]
    set script [file join $temporary worker.xcm]
    set resultFile [file join $temporary result.txt]
    set logFile [file join $temporary worker.log]
    set completed 0
    file mkdir [file join $temporary .xspec cache]
    try {
        save all $snapshot
        set output [open $script w]
        puts $output {query yes}
        puts $output [list lmod wine $::_wine_profile_xspec]
        puts $output "@$snapshot"
        puts $output [list source $::_wine_profile_source]
        puts $output [list _wine_profile_run $number $delta $resultFile [pwd] $trialAngle $fitDelta $toler]
        puts $output {exit}
        close $output
        set childFailed [catch {exec env HOME=$temporary xspec - < $script > $logFile 2>@1} message]
        if {![file exists $resultFile]} {
            _wine_profile_fail "*** WINE likelihood profiling error\n    Profile for p$number produced no result.\n    Inspect $logFile."
        }
        set input [open $resultFile r]
        set result [read $input]
        close $input
        if {[lindex $result 0] ne "ok"} {
            set code [lindex $result 2]
            if {[lrange $code 0 1] eq {WINE BETTER_FIT}} {
                lassign [lrange $code 2 end] angle stat gain threshold name
                _wine_profile_fail [format \
                    "*** WINE likelihood profiling error\n    Profile for p%d stopped: a better fit exceeds toler %.2f.\n    Inclination %.2f deg; statistic %.2f (improvement %.2f).\n    Saved %s in %s. Inspect it or increase toler." \
                    $number $threshold $angle $stat $gain $name [pwd]]
            }
            _wine_profile_fail "*** WINE likelihood profiling error\n    Profile for p$number failed: [lindex $result 1]\n    Diagnostics: $temporary; inspect worker.log or retry the angle with -trial."
        }
        if {$childFailed} {
            _wine_profile_fail "*** WINE likelihood profiling error\n    Profile for p$number stopped in child XSPEC.\n    Inspect $logFile."
        }
        if {$trialAngle ne ""} {
            set ::xspec_tclout [lindex $result 1]
            puts [format "    Trial (p%d): incl %.2f deg, statistic %.2f" \
                $number $trialAngle $::xspec_tclout]
            set completed 1
            return $::xspec_tclout
        }
        lassign [lrange $result 1 end] center bestStat lower upper points
        set ::wine_profile_points $points
        set ::wine_profile_sides [list $lower $upper]
        _wine_profile_present $number $delta $center $lower $upper
        set completed 1
        return $::xspec_tclout
    } finally {
        # Keep the XSPEC transcript and snapshot for diagnosis on any failure.
        if {$completed} {file delete -force $temporary}
    }
}
