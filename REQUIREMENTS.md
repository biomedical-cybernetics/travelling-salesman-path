# Requirements

This software requires [Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) to execute calculations. Follow the instructions below to install this dependency in your operating system.

## Windows

### Cygwin

1. Download and install Cygwin (**32 bit version!**) from [https://www.cygwin.com/](https://www.cygwin.com/).
2. When installing chose a path of your own preference; however, we recommend to make the installation at `C:/cygwin32/bin/`. If you chose a different path, please update the [LocalSettings.m](./LocalSettings.m) file accordingly.

*Note: It is not mandatory to add Cygwin into your environment variables for running concorde.*

### Concorde

1. Go to [https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm](https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm).
2. Download *Concorde for Windows/Cygwin* and *Linkern for Windows/Cygwin*.
3. Extract the executables (they are compressed!).
4. Place the executables in your Cygwin's bin folder. For instance, copy them into `C:/cygwin32/bin/` if you followed the suggested installation in the previous section.
5. Open your *Windows PowerShell* and run the following command to confirm that your installation was successful.

```shell
C:/cygwin32/bin/concorde -h
```

You should see the following output in the console:

```shell
Usage: /usr/bin/concorde [-see below-] [dat_file]
   -B    do not branch
   -C #  maximum chunk size in localcuts (default 16)
   -d    use dfs branching instead of bfs
   -D f  edgegen file for initial edge set
   -e f  initial edge file
   -E f  full edge file (must contain initial edge set)
   -f    write optimal tour as edge file (default is tour file)
   -F f  read extra cuts from file
   -g h  be a grunt for boss h
   -h    be a boss for the branching
   -i    just solve the blossom polytope
   -I    just solve the subtour polytope
   -J #  number of tentative branches
   -k #  number of nodes for random problem
   -K h  use cut server h
   -M f  master file
   -m    use multiple passes of cutting loop
   -n s  problem location (just a name or host:name, not a file name)
   -o f  output file name (for optimal tour)
   -P f  cutpool file
   -q    do not cut the root lp
   -r #  use #x# grid for random points, no dups if #<0
   -R f  restart file
   -s #  random seed
   -S f  problem file
   -t f  tour file (in node node node format)
   -u v  initial upperbound
   -U    do not permit branching on subtour inequalities
   -v    verbose (turn on lots of messages)
   -V    just run fast cuts
   -w    just subtours and trivial blossoms
   -x    delete files on completion (sav pul mas)
   -X f  write the last root fractional solution to f
   -y    use simple cutting and branching in DFS
   -z #  dump the #-lowest reduced cost edges to file xxx.rcn
   -N #  norm (must specify if dat file is not a TSPLIB file)
         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,
         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL
         17=GEOM, 18=JOHNSON
```

If you cannot see this output, then your installation it was not correctly done. Thus, check the steps carefully and try again.

## Linux

1. Go to [https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm](https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm).
2. Download the *Concorde* and *Linkern* unix executables.
3. Create a new `./bin` folder in the root directory of this project. If you chose a different path, please update the [LocalSettings.m](./LocalSettings.m) file accordingly.
4. Place the downloaded executables in the `./bin` folder your created in the previous step.
5. Open your console and run the following command to confirm that your installation was successful.

```shell
./bin/concorde -h
```

You should see the following output in the console:

```shell
Usage: /usr/bin/concorde [-see below-] [dat_file]
   -B    do not branch
   -C #  maximum chunk size in localcuts (default 16)
   -d    use dfs branching instead of bfs
   -D f  edgegen file for initial edge set
   -e f  initial edge file
   -E f  full edge file (must contain initial edge set)
   -f    write optimal tour as edge file (default is tour file)
   -F f  read extra cuts from file
   -g h  be a grunt for boss h
   -h    be a boss for the branching
   -i    just solve the blossom polytope
   -I    just solve the subtour polytope
   -J #  number of tentative branches
   -k #  number of nodes for random problem
   -K h  use cut server h
   -M f  master file
   -m    use multiple passes of cutting loop
   -n s  problem location (just a name or host:name, not a file name)
   -o f  output file name (for optimal tour)
   -P f  cutpool file
   -q    do not cut the root lp
   -r #  use #x# grid for random points, no dups if #<0
   -R f  restart file
   -s #  random seed
   -S f  problem file
   -t f  tour file (in node node node format)
   -u v  initial upperbound
   -U    do not permit branching on subtour inequalities
   -v    verbose (turn on lots of messages)
   -V    just run fast cuts
   -w    just subtours and trivial blossoms
   -x    delete files on completion (sav pul mas)
   -X f  write the last root fractional solution to f
   -y    use simple cutting and branching in DFS
   -z #  dump the #-lowest reduced cost edges to file xxx.rcn
   -N #  norm (must specify if dat file is not a TSPLIB file)
         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,
         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL
         17=GEOM, 18=JOHNSON
```

If you cannot see this output, then check if the downloaded executables have permissions for execution. You can enable them as follows:

```shell
chmod +x ./bin/concorde
chmmd +x ./bin/linkern
```

Then try to execute concorde again.
