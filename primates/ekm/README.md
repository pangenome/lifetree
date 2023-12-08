Dear diary -
---

I am running analyses on the CBCB Nexus Cluster, using nodes with 32 AMD EPYC-7313 cores and 2TB of RAM; see [here](https://wiki.umiacs.umd.edu/umiacs/index.php/Nexus/CBCB).

Interactive jobs
---
```
alias swork="srun --pty --ntasks=1 --mem=4GB --constraint=EPYC-7313 --qos=highmem --partition=cbcb --account=cbcb --time 08:00:00 bash"
swork
```

Installation of seqwish
---

