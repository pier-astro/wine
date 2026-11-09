proc _redwine_profile_check_improvement {angle stat bestStat toler} {
    set gain [expr {$bestStat - $stat}]
    if {$gain <= $toler} {return}
    set index 1
    while {1} {
        set name [expr {$index == 1 ? "better-fit.xcm" : "better-fit-$index.xcm"}]
        set candidate [file join $::_redwine_profile_output_dir $name]
        if {![file exists $candidate]} {break}
        incr index
    }
    if {[catch {save all $candidate} message]} {
        _wind_fail "better fit found, but could not save $candidate: $message"
    }
    return -code error -errorcode [list RED-WINE BETTER_FIT $angle $stat $gain $toler $name] \
        "a better fit was found"
}

proc _redwine_profile_side {ctx center bestStat delta toler direction pointsVar trialCommand} {
    upvar 1 $pointsVar points
    set previous $center
    set previousStat $bestStat
    set step 1.0
    lassign [dict get $ctx profile_bounds] low high
    set edge [expr {$direction < 0 ? $low : $high}]
    if {($direction < 0 && $center <= $edge) || ($direction > 0 && $center >= $edge)} {return [list boundary $center]}
    for {set trial 0} {$trial < 14} {incr trial} {
        set angle [expr {$direction < 0 ? max($edge, $center - $step) : min($edge, $center + $step)}]
        if {abs($angle - $previous) < 1.0e-8} {return [list boundary $previous]}
        set stat [{*}$trialCommand $ctx $angle]
        lappend points [list $angle $stat]
        set excess [expr {$stat - $bestStat}]
        _redwine_profile_check_improvement $angle $stat $bestStat $toler
        if {$excess >= $delta} {
            set inside $previous
            set outside $angle
            set insideStat $previousStat
            set outsideStat $stat
            for {set refine 0} {$refine < 22} {incr refine} {
                set middle [expr {($inside + $outside) / 2.0}]
                set stat [{*}$trialCommand $ctx $middle]
                lappend points [list $middle $stat]
                set excess [expr {$stat - $bestStat}]
                _redwine_profile_check_improvement $middle $stat $bestStat $toler
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
            _wind_fail "profile search did not converge on the [expr {$direction < 0 ? "lower" : "upper"}] side"
        }
        if {$angle == $edge} {return [list boundary $edge]}
        set previous $angle
        set previousStat $stat
        set step [expr {$step * 2.0}]
    }
    _wind_fail "could not bracket the profile limit"
}

