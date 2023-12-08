Dear diary -
---

I am running analyses on the CBCB Nexus Cluster, using nodes with 32 AMD EPYC-7313 cores and 2TB of RAM; see [here](https://wiki.umiacs.umd.edu/umiacs/index.php/Nexus/CBCB).

Interactive jobs
---
```
alias swork="srun --pty --ntasks=1 --mem=4GB --constraint=EPYC-7313 --qos=highmem --partition=cbcb --account=cbcb --time 08:00:00 bash"
```

Installation of seqwish
---
**1.** Start interactive job because need to do build on EPYC-7313.
```
swork
uname -n
```
and check you are on a CBCB node (got `cbcb00.umiacs.umd.edu`).

**2.** Switch gcc version.
```
module load gcc/11.2.0
gcc -v
```
gives
```
Using built-in specs.
COLLECT_GCC=gcc
COLLECT_LTO_WRAPPER=/opt/local/stow/gcc-11.2.0/libexec/gcc/x86_64-pc-linux-gnu/11.2.0/lto-wrapper
Target: x86_64-pc-linux-gnu
Configured with: ./configure --prefix=/opt/local/stow/gcc-11.2.0 --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-languages=c,c++,objc,obj-c++,fortran --disable-dssi --enable-plugin --with-cpu=generic --with-mpfr --with-gmp --with-mpc --with-ppl --with-cloog --disable-multilib
Thread model: posix
Supported LTO compression algorithms: zlib
gcc version 11.2.0 (GCC)
```
Check you can find `libdl` by typing `whereis libdl`, which should give 
```
libdl: /usr/lib64/libdl.so /usr/lib64/libdl.a
```
Check you can find `zlib` by typing `whereis zlip`, which should give
```
zlib: /usr/include/zlib.h /usr/share/man/man3/zlib.3.gz
```
(Note you sometimes need to load `zlib` with command `module avail zlib`).

**3.** Make path to software directory
```
SWDIR=/fs/cbcb-lab/ekmolloy/ekmolloy/lifetree/primates/ekm/software
```

**4.** Install [jemalloc](https://jemalloc.net). Note last commit was [e4817c8](https://github.com/jemalloc/jemalloc/commit/e4817c8d89a2a413e835c4adeab5c5c4412f9235).
```
cd $SWDIR
git clone https://github.com/jemalloc/jemalloc.git
mkdir jemalloc-install
cd jemalloc
./autogen.sh
./configure --prefix="$SWDIR/jemalloc-install"
make && make install
```
and then check static lib has been created here: `$SWDIR/jemalloc-install/lib`.

**5.** Install [seqwish](https://github.com/ekg/seqwish.git). Note that last commit was [f44b402](https://github.com/ekg/seqwish/commit/f44b402f0c2e02988d431d9b2e5eba9727cf93a9).
```
cd  $SWDIR
git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
cmake -DCMAKE_INSTALL_PREFIX=/ -H. -Bbuild && cmake --build build -- -j 3
```
This will fail so you need to `rm -rf build`, update [CMakeList](), and try again.
