#!/usr/bin/env python3
"""Live Tcl geometry/scan check with an isolated, file-backed toy spectrum.

Run with initialized PyXSPEC and XSPEC on PATH. Uses a unit-area diagonal
response so native session snapshots can reload it (unlike a dummy response).
No observational data, Poisson significance, or global recovery test.
"""
from pathlib import Path
import os
import subprocess
import tempfile

import numpy as np
from astropy.io import fits
from xspec import AllData, AllModels, FakeitSettings, Fit, Model, Xset


def main():
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="wind-procedures-") as temporary:
        work = Path(temporary)
        edges = np.linspace(.5, 12., 121).astype("f4")
        channels = np.arange(1, 121, dtype="i2")
        ones = np.ones(120, dtype="i2")
        matrix = fits.BinTableHDU.from_columns([
            fits.Column(name="ENERG_LO", format="E", unit="keV", array=edges[:-1]),
            fits.Column(name="ENERG_HI", format="E", unit="keV", array=edges[1:]),
            fits.Column(name="N_GRP", format="I", array=ones),
            fits.Column(name="F_CHAN", format="I", array=channels),
            fits.Column(name="N_CHAN", format="I", array=ones),
            fits.Column(name="MATRIX", format="E", array=np.ones(120, dtype="f4")),
        ], name="MATRIX")
        ebounds = fits.BinTableHDU.from_columns([
            fits.Column(name="CHANNEL", format="I", array=channels),
            fits.Column(name="E_MIN", format="E", unit="keV", array=edges[:-1]),
            fits.Column(name="E_MAX", format="E", unit="keV", array=edges[1:]),
        ], name="EBOUNDS")
        for hdu, kind in ((matrix, "RSP_MATRIX"), (ebounds, "EBOUNDS")):
            hdu.header.update(HDUCLASS="OGIP", HDUCLAS1="RESPONSE", HDUCLAS2=kind,
                              HDUVERS="1.3.0", DETCHANS=120, CHANTYPE="PI",
                              TELESCOP="TEST", INSTRUME="TEST", FILTER="NONE")
        matrix.header.update(HDUCLAS3="FULL", TLMIN4=1)
        ebounds.header["TLMIN1"] = 1
        response = work / "unit.rmf"
        fits.HDUList([fits.PrimaryHDU(), matrix, ebounds]).writeto(response)

        Xset.chatter = 0
        Xset.allowPrompting = False
        AllModels.lmod("wine", str(root / "xspec"))
        AllData.dummyrsp(.5, 12., 120, "lin", .5, 11.5/120)
        model = Model("powerlaw+windemlos*gaussian+windemlos*gaussian")
        model.setPars(2., .01, .25, .4396926207859083, .242050121690, 35.,
                      6.7, .025, .005, .25, .4396926207859083, .242050121690,
                      35., 6.7, .025, .0025)
        for slave, master in ((10, 3), (11, 4), (12, 5), (13, 6)):
            model(slave).link = f"p{master}"
        AllModels.setEnergies("0.1 30 2000 log")
        AllData.fakeit(1, FakeitSettings(response=str(response), exposure=100000,
                                       fileName=str(work / "toy.fak")), applyStats=False)
        AllData(1).response = str(response)
        Fit.statMethod = "cstat"
        Fit.nIterations = 40
        Fit.criticalDelta = 1e-6
        for index in (1, 3, 7, 8, 9, 14, 15, 16):
            model(index).frozen = True
        snapshot = work / "reference.xcm"
        Xset.save(str(snapshot))
        with snapshot.open("a") as stream:
            stream.write("\nenergies 0.1 30 2000 log\n")
        # Fixed off-wind geometry permits an exact native view->incl transform.
        AllModels.clear()
        off = Model("powerlaw+windemoff*gaussian")
        hole = (1-np.cos(np.deg2rad(20.)))/(1-.4396926207859083)
        off.setPars(2., .01, .25, .4396926207859083, hole, .5, 0,
                    6.7, .025, .0075)
        for index in (1, 3, 4, 5, 7, 8, 9, 10):
            off(index).frozen = True
        AllData.fakeit(1, FakeitSettings(response=str(response), exposure=100000,
                                       fileName=str(work / "off.fak")), applyStats=False)
        AllData(1).response = str(response)
        Xset.save(str(work / "off.xcm"))
        with (work / "off.xcm").open("a") as stream:
            stream.write("\nenergies 0.1 30 2000 log\n")
        AllModels.clear()
        Xset.restore(str(snapshot))
        # A named model shared by two data groups exercises qualified IDs and
        # selection of an identity-linked slave component in group 2.
        AllData(f"2:2 {work / 'toy.fak'}")
        AllData(2).response = str(response)
        AllModels.clear()
        named = Model("powerlaw+windemlos*gaussian+windemlos*gaussian", "named", 1)
        named.setPars(2., .01, .25, .4396926207859083, .242050121690, 35.,
                      6.7, .025, .005, .25, .4396926207859083, .242050121690,
                      35., 6.7, .025, .0025)
        for slave, master in ((10, 3), (11, 4), (12, 5), (13, 6)):
            named(slave).link = f"named:p{master}"
        for index in (1, 3, 7, 8, 9, 14, 15, 16):
            named(index).frozen = True
        Xset.save(str(work / "named.xcm"))
        with (work / "named.xcm").open("a") as stream:
            stream.write("\nenergies 0.1 30 2000 log\n")

        # No test hooks in the public command: failure injection affects only
        # this disposable native process and its internal worker.
        script = r'''
query no
lmod wine @ROOT@/xspec
@@WORK@/reference.xcm
source @ROOT@/procedures/wind_geometry.tcl
if {[catch {
    if {[windgeom 4] ne ""} {_wind_fail "Structured result echoed instead of stored"}
    set noCovariance $::xspec_tclout
    if {[dict get $noCovariance uncertainty_status] ne "unavailable"} {_wind_fail "Unexpected covariance before fitting"}
    fit
    save all @WORK@/before.xcm
    set beforeGrid [_wind_out energies 1]
    set beforeStat [_wind_out stat]
    windgeom 4
    set geometry $::xspec_tclout
    if {abs([dict get $geometry tin]-20)>1e-3} {_wind_fail "Wrong linked component"}
    source @FITEXPORT@
    set sigma [coverror {$p(4)}]
    if {abs($sigma-[dict get $geometry uncertainty Cf])>1e-6} {_wind_fail "Covariance Cf mismatch"}
    proc ::tcl::mathfunc::windtin {cf skew incl} {
        return [dict get [_wind_angles [dict create kind los] [dict create Cf $cf skew $skew incl $incl]] tin]
    }
    set sigma [coverror {windtin($p(4),$p(5),$p(6))}]
    if {abs($sigma-[dict get $geometry uncertainty tin])>1e-4*max(1.0,$sigma)} {_wind_fail "Covariance angular cross-term mismatch"}
    foreach quantity {tin tout incl} {
        if {[winderror 4 $quantity delta 1.0 toler 1.0] ne ""} {_wind_fail "Error result echoed"}
        set limits $::xspec_tclout
        if {[dict get $limits lower]>=[dict get $limits center] ||
            [dict get $limits upper]<=[dict get $limits center]} {_wind_fail "Limits do not bracket the reference"}
        if {[dict get $limits flag_format] eq "wind"} {
            foreach side [dict get $limits sides] {
                if {[lindex $side 0] ne "ok"} {_wind_fail "Unexpected boundary/unstable toy limit"}
            }
        }
    }
    if {[windsteppar 4 tin 0 20 1 tout 60 90 1] ne ""} {_wind_fail "Scan result echoed"}
    set scan $::xspec_tclout
    foreach row [dict get $scan rows] {
        if {[dict get $row status] ne "ok"} {_wind_fail "Permitted fit failed"}
    }
    windsteppar 2 tin 0 60 1 tout 20 60 1
    set scan $::xspec_tclout
    save all @WORK@/after.xcm
    if {$beforeGrid ne [_wind_out energies 1] || $beforeStat ne [_wind_out stat]} {
        _wind_fail "Caller energies/statistic changed"
    }
    source @EXPORT@
    savesteppar result $scan outfile @WORK@/scan.fits clobber yes
    freeze 4
    winderror 4 tin delta 1.0 toler 1.0
    set fixedCover $::xspec_tclout
    if {[dict get $fixedCover status] ne "ok"} {_wind_fail "Fixed-Cf profile failed"}
    freeze 5
    winderror 4 tout delta 1.0 toler 1.0
    set fixedShape $::xspec_tclout
    if {[dict get $fixedShape mapping] ne "monotonic-fixed-geometry"} {_wind_fail "Missing fixed-geometry transformation"}
    thaw 5
    thaw 4
    steppar best 4 0.4 0.45 1 6 30 40 1
    savesteppar 4 6 outfile @WORK@/native.fits clobber yes
    # Exercise the real worker's Tcl error path, preserving excluded cells.
    rename fit _real_fit
    proc fit {args} {return -code error "simulated fit failure"}
    _wind_scan_worker 2 [_wind_options] [dict get $scan axes] @WORK@/reference.xcm @WORK@
    set handle [open @WORK@/result.tcl r]
    set failedScan [read $handle]
    close $handle
    set failures 0
    foreach row [dict get $failedScan rows] {
        if {[dict get $row status] eq "failed"} {
            if {[dict get $row message] ne "simulated fit failure"} {_wind_fail "Lost fit error"}
            incr failures
        }
    }
    if {$failures != 3} {_wind_fail "Wrong failed/excluded accounting"}
    rename fit {}
    rename _real_fit fit
    # A narrowed native inclination range may touch the requested shell at
    # exactly one point. Freeze that inclination instead of setting equal bounds.
    newpar 6 20 0.01 0 0 20 20
    set touchingContext [_wind_context 2 [_wind_options]]
    _wind_constrain $touchingContext 20 60
    set touched [_wind_angles $touchingContext]
    if {abs([dict get $touched incl]-20)>1e-6} {_wind_fail "Incorrect one-point inclination range"}
    model clear
    @@WORK@/named.xcm
    fit
    windgeom 4 model named group 2
    set geometry $::xspec_tclout
    if {abs([dict get $geometry tin]-20)>1e-3} {_wind_fail "Wrong named/group geometry"}
    if {[dict get $geometry uncertainty_status] ne "local"} {_wind_fail "Missing named/group covariance"}
    winderror 4 incl delta 1.0 model named group 2
    set namedLimits $::xspec_tclout
    if {[dict get $namedLimits lower]>=35 || [dict get $namedLimits upper]<=35} {_wind_fail "Invalid named/group inclination limits"}
    windsteppar 4 tin 0 20 1 tout 60 90 1 model named group 2
    set namedScan $::xspec_tclout
    foreach row [dict get $namedScan rows] {
        if {[dict get $row status] ne "ok"} {_wind_fail "Named/group scan failed"}
    }
    model clear
    @@WORK@/off.xcm
    fit
    windgeom 2
    set offGeometry $::xspec_tclout
    winderror 2 incl delta 1.0
    set offLimits $::xspec_tclout
    foreach side {lower upper} raw [dict get $offLimits native_limits] {
        set expected [expr {$raw == 0 ? 0 : $raw*[dict get $offGeometry tin]}]
        if {abs([dict get $offLimits $side]-$expected)>1e-6} {_wind_fail "Incorrect off-wind physical limit"}
    }
    thaw 4 5
    if {![catch {winderror 2 incl delta 1.0} rejected] || ![string match {*require fixed Cf and hole*} $rejected]} {
        _wind_fail "Unconstrained off-wind inclination must not be relabelled from view"
    }
    foreach quantity {tin tout} target {25.0 65.0} {
        @@WORK@/off.xcm
        thaw 4 5
        set context [_wind_context 2 [_wind_options]]
        _wind_single_constraint $context $quantity $target
        set constrained [_wind_angles $context]
        if {abs([dict get $constrained $quantity]-$target)>1e-5} {_wind_fail "Incorrect native hole/Cf angular link"}
    }
    untie 4
    newpar 4 1.0
    freeze 4
    if {![catch {winderror 2 tin} rejected] || ![string match {*fixed*full-coverage*} $rejected]} {
        _wind_fail "Frozen full-coverage opening angles must be reported as fixed"
    }
} message options]} {
    puts "TEST FAILED: $message"
    puts [dict get $options -errorinfo]
} else {puts "WIND PROCEDURES PASSED"}
exit
'''.replace("@ROOT@", str(root)).replace("@WORK@", str(work)).replace(
            "@EXPORT@", str(root.parent / "xspeclab/procedures/export_steppar.tcl")).replace(
            "@FITEXPORT@", str(root.parent / "xspeclab/procedures/export_fit.tcl"))
        # The native dummy response is not involved in any reloaded session.
        environment = os.environ.copy()
        (work / "pfiles").mkdir()
        environment["PFILES"] = f"{work / 'pfiles'};{environment['HEADAS']}/syspfiles"
        environment["HEADASNOQUERY"] = "1"
        native_script = work / "check.tcl"
        native_script.write_text(script)
        completed = subprocess.run(["xspec", "-"], input=f"source {native_script}\nexit\n", text=True,
                                   capture_output=True, cwd=work, env=environment,
                                   timeout=180)
        log = completed.stdout + completed.stderr
        assert "\nWIND PROCEDURES PASSED\n" in log, log[-10000:]
        assert (work / "after.xcm").exists(), log[-18000:]
        assert (work / "before.xcm").read_text() == (work / "after.xcm").read_text()
        with fits.open(work / "scan.fits") as hdus:
            table = hdus[1]
            assert table.header["PNAME1"] == "tin"
            assert "PIDX1" not in table.header
            assert np.sum(table.data["STATUS"] == "excluded") == 1
            assert np.sum(np.isfinite(table.data["STAT"])) == 3
            assert np.isnan(table.data["STAT"][table.data["STATUS"] == "excluded"]).all()
        with fits.open(work / "native.fits") as hdus:
            assert hdus[1].header['PNAME1'] == 'Cf'
            assert hdus[1].header['PIDX1'] == 4
            assert len(hdus[1].data) == 4
        print("Native angular profiles, covariance versus coverror, linked scans,")
        print("named/group contexts, boundaries, session preservation,")
        print("failed/excluded statuses and FITS export passed")


if __name__ == "__main__":
    main()
