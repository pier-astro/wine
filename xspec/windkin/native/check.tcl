# Run in a fresh XSPEC session after lmod redwine and lmod windkin.
dummyrsp 0.2 20 1000 log
foreach expression {wuline wgline zwgline wuconv*gaussian wgconv*gaussian wglinelos_Cf zwglinelos_Cf wgconvlos_Cf*gaussian zwline wsectline zwsectline wsect*gaussian} {
    model $expression & /*
    flux 0.2 20
    tclout flux 1
    set photons [lindex $xspec_tclout 3]
    if {!($photons > 0 && $photons < 100)} {
        return -code error "Invalid model flux for $expression: $xspec_tclout"
    }
}
foreach {expression redshift_index} {zwline 6 zwgline 8 zwsectline 9} {
    dummyrsp 0.2 20 1000 log
    model $expression & /*
    flux 0.2 20
    tclout flux 1
    set rest [lindex $xspec_tclout 3]
    newpar $redshift_index 0.184
    # Match the source-frame bins used for the reference line.
    dummyrsp [expr {0.2/1.184}] [expr {20.0/1.184}] 1000 log
    flux 0.1 20
    tclout flux 1
    set shifted [lindex $xspec_tclout 3]
    if {abs(1.184*$shifted/$rest - 1.0) > 1e-4} {
        return -code error "Unexpected integrated redshift scaling for $expression: rest=$rest shifted=$shifted"
    }
}
# Reflecting inclination and the signed PDF gives the same profile.
model zwgline & /*
newpar 1 -0.02 0.001 -0.99 -0.99 0.99 0.99
newpar 2 0.08
flux 0.1 20
tclout flux 1
set inward [lindex $xspec_tclout 3]
newpar 1 0.02
newpar 5 135
flux 0.1 20
tclout flux 1
set opposite [lindex $xspec_tclout 3]
if {abs($inward/$opposite - 1.0) > 1e-7} {
    return -code error "Negative vmean does not match reversed outflow"
}
model wsectline & /*
newpar 1 0
newpar 6 360
flux 0.1 20
tclout flux 1
set full [lindex $xspec_tclout 3]
newpar 6 180
flux 0.1 20
tclout flux 1
set half [lindex $xspec_tclout 3]
if {abs(2.0*$half/$full - 1.0) > 1e-4} {
    return -code error "Static half-sector flux is $half; expected half of $full"
}
puts "windkin native checks passed"
