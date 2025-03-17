Before running the program, load all the necessary modules on Teach cluster.
```
source teachsetup
```

Compile the program. This command will output an executable `stft`.
```
make
```

`stft` takes two arguments, `num-samples` and `window-size`, in this order. It also requires options for the algorithm `-a` and the parallel scheme `-p`. For detailed usage of `stft`, run `stft -h`. 

Performance analysis done on Teach cluster are in `/presentation-analysis` and `/report-analysis`. The jobscripts and plot scripts used for analysis are also provided.