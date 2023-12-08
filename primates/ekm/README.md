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
1. Start interactive job because need to do build on EPYC-7313.
```
swork
uname -n
```
and check you are on a CBCB node (got `cbcb00.umiacs.umd.edu`).

2. Switch gcc version.
```
module load gcc/11.2.0
gcc -v
```
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

3. 
```
```
