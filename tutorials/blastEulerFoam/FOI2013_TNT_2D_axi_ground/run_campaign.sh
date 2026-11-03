#!/bin/bash
#------------------------------------------------------------------------------
# run_campaign.sh: several FOI 2013 shots with run_foi2013.sh, at most -jobs
# at a time (each 1D stage on one core, each 2D stage on -np cores).
#
#   run_campaign.sh [-jobs 2] [-np 6] [-mixingLevel 3] [-dir foi2013_runs]
#                   [shot:closure ...]
#
# Default list, in order of priority: HOB 2 m (4 TNT, 5 TNT/Al), 0.5 m
# (2 and 9 TNT repeats, 7 TNT/Al), 1.0 m (3, 6), 0.15 m (1, 8), laminar
# first, then PaSR. One log per run: <dir>/run_shot<N>_<closure>.log.
# Run detached:
#   setsid nohup ./run_campaign.sh -jobs 4 -np 8 > campaign.log 2>&1 &
#------------------------------------------------------------------------------
HERE=$(cd "${0%/*}" && pwd -P)
jobs=2 np=6 mixingLevel=3 dir=foi2013_runs runs=
while [ $# -gt 0 ]; do
    case "$1" in
        -jobs) jobs=$2; shift ;;
        -np) np=$2; shift ;;
        -mixingLevel) mixingLevel=$2; shift ;;
        -dir) dir=$2; shift ;;
        -h|-help) sed -n '3,14p' "$0"; exit 0 ;;
        *:*) runs="$runs $1" ;;
        *) echo "unknown argument $1 (runs are shot:closure)" >&2; exit 1 ;;
    esac
    shift
done
if [ -z "$runs" ]; then
    for c in laminar pasrd; do
        for s in 4 5 2 9 7 3 6 1 8; do runs="$runs $s:$c"; done
    done
fi
mkdir -p "$dir"
dir=$(cd "$dir" && pwd -P)
echo "$(date): $jobs at a time, 2D on $np processes, mixing level $mixingLevel:$runs"
printf '%s\n' $runs | xargs -P "$jobs" -I{} bash -c '
    s=${1%%:*}; c=${1##*:}
    if "$0/run_foi2013.sh" -shot "$s" -closure "$c" -np "$2" -mixingLevel "$3" -dir "$4" \
        > "$4/run_shot${s}_${c}.log" 2>&1
    then echo "$(date): shot $s $c done"
    else echo "$(date): shot $s $c FAILED, see $4/run_shot${s}_${c}.log"
    fi' "$HERE" {} "$np" "$mixingLevel" "$dir"
echo "$(date): campaign finished"

#------------------------------------------------------------------------------
