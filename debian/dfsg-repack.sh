#!/bin/bash
# Borrowed from Raphael Geissert's Debian PHP repack script.

set -e

if [ ! -f "$3" ] && [ ! -f "$1" ]; then
    echo "This script must be run via uscan or by manually specifying the tarball" >&2
    exit 1
fi

tarball=

[ -f "$3" ] && tarball="$3"
[ -z "$tarball" -a -f "$1" ] && tarball="$1"

fname="$(basename "$tarball")"
tarball="$(readlink -f "$tarball")"

tdir="$(mktemp -d)"
trap '[ ! -d "$tdir" ] || rm -r "$tdir"' EXIT

xzcat "$tarball" | tar --wildcards --delete '*/mex_interface/*.mat' --delete '*/examples/lib_win64' --delete '*/examples/*_win64.*' --delete '*/*pdf' > "$tdir/${fname/.xz}"
#touch -m -r "$tarball" "$tdir/${fname/.gz}"
xz "$tdir/${fname/.xz}"

mv "$tarball" "$tarball.bkp"
tarballdfsg=$(echo $3|sed 's/\.orig\.tar\.xz/+dfsg.orig.tar.xz/')
echo Writing $tarballdfsg
mv "$tdir/$fname" "$tarballdfsg"
