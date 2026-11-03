#!/bin/bash
#------------------------------------------------------------------------------
# run_foi2013.sh: one FOI June 2013 shot, 1D free-air stage then the 2D
# axisymmetric case with the ground.
#
#   run_foi2013.sh -shot N [-closure laminar|pasrd] [-np 6] [-mixingLevel 3]
#                  [-dir <work dir>] [-only1D] [-t0 <mapping time>]
#                  [-end1D 0.005]
#
# Work directory <dir>/shot<N>_<closure>/ with 1D/ and 2D/ (default dir:
# ./foi2013_runs). Each stage is skipped if it has finished ("End" in its
# solver log; the 1D also at 5 ms), so the script can be rerun after an
# interruption; a stage that stopped part way is started again from
# scratch (its directory is deleted).
#
# 1D (serial): setup_foi2013.py, run to the mapping time t0 (T0_BY_HOB in
# setup_foi2013_2d.py), so that t0 is written exactly, then continued to
# 5 ms (-end1D; the outer probes: the 2D level table needs their arrival
# times).
# 2D (-np processes): setup_foi2013_2d.py and its Allrun (mapping at t0,
# then blastEulerFoam to 18 ms).
#
# Requirements: OpenFOAM-12 and this blastfoam build sourced, Python 3 with
# numpy, the EuDetCode fits in utilities/thermo_fitter/eudet_foi2013 (all
# shots are in the repository). Run detached, e.g.
#   setsid nohup ./run_foi2013.sh -shot 4 -np 6 > run_shot4.log 2>&1 &
#------------------------------------------------------------------------------
set -e
HERE=$(cd "${0%/*}" && pwd -P)
BLAST=${BLAST_DIR:-$(cd "$HERE/../../.." && pwd -P)}
ONED_SETUP="$BLAST/tutorials/blastEulerFoam/FOI2013_TNT_1D_spherical/setup_foi2013.py"
TWOD_SETUP="$HERE/setup_foi2013_2d.py"

shot= closure=laminar np=6 mixingLevel=3 dir=foi2013_runs only1D= t0= end1D=0.005
while [ $# -gt 0 ]; do
    case "$1" in
        -shot) shot=$2; shift ;;
        -closure) closure=$2; shift ;;
        -np) np=$2; shift ;;
        -mixingLevel) mixingLevel=$2; shift ;;
        -dir) dir=$2; shift ;;
        -only1D) only1D=yes ;;
        -t0) t0=$2; shift ;;
        -end1D) end1D=$2; shift ;;
        -h|-help) sed -n '3,25p' "$0"; exit 0 ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
    shift
done
[ -n "$shot" ] || { echo "-shot N required" >&2; exit 1; }
command -v blastEulerFoam >/dev/null || { echo "blastEulerFoam not found: source OpenFOAM-12 and blastfoam" >&2; exit 1; }
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

W=$(mkdir -p "$dir/shot${shot}_${closure}" && cd "$dir/shot${shot}_${closure}" && pwd -P)
if [ -z "$t0" ]; then
    t0=$(cd "$HERE" && python3 -c "
import sys; sys.path.insert(0, '$BLAST/utilities/thermo_fitter')
import foi2013_charges as fc, setup_foi2013_2d as s
print(s.T0_BY_HOB[fc.SHOTS['$shot'][2]])")
fi
echo "shot $shot, closure $closure: work dir $W, mapping time $t0"

finished() { [ -f "$1/log.blastEulerFoam" ] && grep -q "^End" "$1/log.blastEulerFoam"; }
# the 1D log has an End after each leg: also require the final time
finished1D() {
    finished "$1" && [ -d "$1/$t0" ] && grep -a "^Time = " "$1/log.blastEulerFoam" | tail -1 \
        | awk -v e="$end1D" '{exit !($3 >= e*(1 - 1e-6))}'
}

# -- 1D free-air stage
oneD="$W/1D"
if finished1D "$oneD"; then
    echo "1D: finished"
else
    rm -rf "${oneD:?}"
    python3 "$ONED_SETUP" --shot "$shot" --closure "$closure" --out "$oneD" --endTime "$end1D"
    (
        cd "$oneD"
        foamDictionary -entry endTime -set "$t0" system/controlDict >/dev/null
        foamDictionary -entry writeInterval -set "$t0" system/controlDict >/dev/null
        echo "1D: to t0 = $t0 ($(date))"
        ./Allrun > allrun.out 2>&1
        [ -d "$t0" ] || { echo "1D: no time directory $t0" >&2; exit 1; }
        # Continued from t0: <t0>/uniform/time keeps beginTime 0, so this is a
        # restart (the detonation energy is not added again)
        foamDictionary -entry startFrom -set latestTime system/controlDict >/dev/null
        foamDictionary -entry endTime -set "$end1D" system/controlDict >/dev/null
        foamDictionary -entry writeInterval -set 5e-5 system/controlDict >/dev/null
        echo "1D: continued to $end1D s ($(date))"
        runApplication -a blastEulerFoam
        python3 summary.py > summary.txt
    )
    finished1D "$oneD" || { echo "1D did not finish: $oneD/log.blastEulerFoam" >&2; exit 1; }
fi
[ -n "$only1D" ] && exit 0

# -- 2D case with the ground
twoD="$W/2D"
if finished "$twoD"; then
    echo "2D: finished"
    exit 0
fi
rm -rf "${twoD:?}"
python3 "$TWOD_SETUP" --oneD "$oneD" --out "$twoD" --t0 "$t0" --np "$np" \
    --mixingLevel "$mixingLevel"
echo "2D: mapping and run on $np processes ($(date))"
(cd "$twoD" && ./Allrun > allrun.out 2>&1)
finished "$twoD" && echo "2D: finished ($(date))" || { echo "2D did not finish: $twoD" >&2; exit 1; }

#------------------------------------------------------------------------------
