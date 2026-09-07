#!/usr/bin/env bash
# Cross-build ballistics-engine for riscv64 on bosgame (amd64).
#
#   linux   -> riscv64gc-unknown-linux-gnu    tier 2: prebuilt std, GNU cross gcc
#   freebsd -> riscv64gc-unknown-freebsd      tier 3: -Z build-std + harvested sysroot
#   netbsd  -> riscv64gc-unknown-netbsd       tier 3: ditto
#   openbsd -> riscv64gc-unknown-openbsd      tier 3: ditto
#
# The BSD targets have no prebuilt std, so they need nightly and a sysroot taken
# from the matching verification VM. Linux deliberately uses gnu, not musl: the
# shipped artifact is dynamically linked against glibc, and a musl build fails
# anyway when printpdf links a shared object with no libc.so in the sysroot.
set -euo pipefail
OS="${1:?usage: cross-riscv64.sh <linux|freebsd|netbsd|openbsd>}"
SRC="${SRC:-/work/ballistics-engine}"
SYSROOTS="${SYSROOTS:-/work/sysroots}"

case "$OS" in
  linux)   TRIPLE=riscv64gc-unknown-linux-gnu ;;
  freebsd) TRIPLE=riscv64gc-unknown-freebsd   ;;
  netbsd)  TRIPLE=riscv64gc-unknown-netbsd    ;;
  openbsd) TRIPLE=riscv64gc-unknown-openbsd   ;;
  *) echo "unknown os: $OS" >&2; exit 2 ;;
esac
ENV_TRIPLE="${TRIPLE//-/_}"

cd "$SRC"
if [ "$OS" = linux ]; then
    rustup target add "$TRIPLE" >/dev/null
    export CARGO_TARGET_${ENV_TRIPLE^^}_LINKER=riscv64-linux-gnu-gcc
    export CC_${ENV_TRIPLE}=riscv64-linux-gnu-gcc
    export AR_${ENV_TRIPLE}=riscv64-linux-gnu-ar
    BUILD=(cargo build --locked --release --target "$TRIPLE")
else
    SYSROOT="$SYSROOTS/$OS"
    [ -d "$SYSROOT" ] || { echo "missing sysroot: $SYSROOT" >&2; exit 1; }
    # cc-rs compiles ring's C; clang cross-targets natively given a sysroot.
    export CC_${ENV_TRIPLE}=clang
    EXTRA_CFLAGS=""
    if [ "$OS" = netbsd ]; then
        # NetBSD's <sys/common_wchar_limits.h> #errors unless __WCHAR_MIN__ and
        # __WINT_MIN__ are predefined. GCC defines them; clang defines only the
        # __*_MAX__ halves, at every version through 19, so the header rejects a
        # perfectly good clang. Supply the missing halves derived from clang's own
        # maxima for this ABI: __WCHAR_MAX__ 2147483647 (signed 32-bit wchar_t)
        # and __WINT_MAX__ 4294967295U (unsigned 32-bit wint_t).
        EXTRA_CFLAGS="-D__WCHAR_MIN__=(-2147483647-1) -D__WINT_MIN__=0U"
    fi
    export CFLAGS_${ENV_TRIPLE}="--target=riscv64-unknown-${OS} --sysroot=${SYSROOT} ${EXTRA_CFLAGS}"
    export AR_${ENV_TRIPLE}=llvm-ar
    export CARGO_TARGET_${ENV_TRIPLE^^}_LINKER=clang
    # --sysroot alone is not enough on every BSD: clang's FreeBSD driver derives
    # the library search path from it, its NetBSD driver does not, so lld reports
    # "unable to find library -lc" against a sysroot that plainly has one. Pass
    # both lib dirs explicitly -- harmless where the driver already found them.
    export CARGO_TARGET_${ENV_TRIPLE^^}_RUSTFLAGS="-C link-arg=--target=riscv64-unknown-${OS} -C link-arg=--sysroot=${SYSROOT} -C link-arg=-fuse-ld=lld -C link-arg=-L${SYSROOT}/usr/lib -C link-arg=-L${SYSROOT}/lib"
    # Tier 3: no prebuilt std, so compile it from source alongside the crate.
    # Use the toolchain already selected, do NOT say `+nightly`: in a nightly
    # image that resolves to a DIFFERENT, freshly-downloaded nightly without the
    # rust-src component, and build-std then fails on a missing Cargo.lock.
    # Substitute, don't pipe: `rustc -vV | grep -q` dies under `set -o pipefail`
    # because grep -q closes the pipe on first match and rustc takes SIGPIPE, so
    # the pipeline reports failure on the very toolchain it just matched.
    [[ "$(rustc -vV)" == *nightly* ]] || { echo "build-std needs a nightly toolchain" >&2; exit 1; }
    BUILD=(cargo build -Z build-std=std,panic_abort --locked --release --target "$TRIPLE")
fi

echo "==> $OS  ($TRIPLE)"
"${BUILD[@]}"
OUT="${CARGO_TARGET_DIR:-target}/$TRIPLE/release/ballistics"
ls -l "$OUT"; file "$OUT"
