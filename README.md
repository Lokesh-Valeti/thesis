# PRAC: Private Random Access Computations

Ian Goldberg, iang@uwaterloo.ca

PRAC implements three-party secure computation, with a particular focus on computations that require random access to memory.  Parties 0 and 1 are the computational peers, while party 2 is the server.  The server aids the computation, but generally does much less than the two computational peers.

The multi-party computation (MPC) makes use of _resources_, most notably multiplication triples and distributed point functions (DPFs).  These resources can be precomputed; they are independent of the values in the computation being performed, so you only need to know how many of each you'll need.

PRAC has three _modes_:

  - Precomputation mode
    - In this mode, the parties precompute all of the resources they will need in the online phase.  The resources will be stored in files on disk.

  - Online mode
    - In this mode, the parties perform the actual computation, using the resources stored by an earlier precomputation mode.
    - **Note:** using resources during online mode does _not_ currently delete them from disk; this is to facilitate timing tests that run the online mode many times, for example with differing numbers of threads or networking configurations.  This is insecure when used in practice, however; be sure to never use a precomputed resource for more than one computation, or private information may be leaked!

  - Online-only mode
    - In this mode, there is no precomputation; all of the required resources are computed on the fly, as needed.  At the end of the run, the number of each resource so computed will be output, so that you can feed that into a precomputation mode for next time.

PRAC supports two kinds of threading: communication threads and local processing threads.

  - Communication threads allow multiple independent subproblems to be performed at the same time, where each subproblem requires multi-party interaction.
  - Local processing threads allow individual parties to perform large local computations (for example, computing large local dot products) in a multithreaded manner, without interacting with other parties.

Currently, all of the interesting multithreading (except for some unit tests) are local processing threads.

## Non-docker instructions:

Build with `make`.  The build has been tested on Ubuntu 20.04 and Ubuntu
22.04.  You'll need libbsd-dev and libboost-all-dev.

You'll need three terminals with the built binary, either on the same machine, or different machines.  If it's on the same machine, it's OK for them all to be in the very same directory; the saved resources for the different parties won't clobber each other.

In each terminal, run <code>./prac _opts_ _player\_num_ _player\_addresses_ _args_</code>

  - <code>_player\_num_</code> is `0`, `1`, or `2`
  - <code>_player\_addresses_</code> is:
    - for player 0: nothing (just leave it out)
    - for player 1: player 0's hostname
    - for player 2: player 0's hostname followed by player 1's hostname
    - (you can use IP addresses instead of hostnames if you want)

For preprocessing mode:

  - <code>_opts_</code> should be `-p` to indicate preprocessing mode
  - You can also add `-a` to indicate that the resources being computed should be _appended_ to existing resource files, rather than replacing them.  This is useful if you want to create for example a large number of large DPFs.  Each DPF creation consumes a lot of memory, and all the DPFs are computed simultaneously in order to save network round trips.  To compute them in smaller batches if your memory is limited, run `./prac` several times, each with the `-a` option.
  - The `-e` option means to store precomputed DPFs in expanded form.  This can take a lot of disk space for large DPFs, and depending on your disk speed and CPU capabilities and number of local processing threads, it can actually be faster to recompute the expansion than to read it from disk.
  - Use the <code>-t _num_</code> option to enable <code>_num_</code> _communication_ threads during preprocessing mode.  As above, this is probably not what you want.

Then the <code>_args_</code> specify what resources to compute.  You actually only need to provide the args to player 2 (the server) in preprocessing mode; it will inform the other players what is being computed.  The other players will ignore their arguments in preprocessing mode.

The args in preprocessing mode are each of the form <code>_type_:_num_</code> to indicate to create <code>_num_</code> resources of the given type.  The available types include:

  - `t`: a multiplication triple
  - `h`: a multiplication half-triple
  - <code>r*d*</code>: A DPF of depth <code>_d_</code> for random accesses to memory (RDPF).  A DPF of depth _d_ can be used to process 2<sup>_d_</sup> memory locations.
  - `c`: a DPF for comparisons (CDPF)

If you do have multiple communication threads, you can optionally specify different sets of resources to create in each one by prefixing communication thread _i_'s list with <code>T_i_</code> (_i_ starts at 0).  If the args do not start with such an indication, each communication thread will use a copy of the provided list.

You can also include the pseudo-resource argument <code>p:_num_</code> to indicate that you want to use <code>_num_</code> local processing threads for this communication thread.

Example:  In three terminals on the same host, run:

  - `./prac -p 0`
  - `./prac -p 1 localhost`
  - `./prac -p 2 localhost localhost t:100 r6:10 r20:5 c:50 p:8`

to create 100 multiplication triples, 10 RDPFs of depth 6, 5 RDPFs of depth 20, and 50 CDPFs, using 8 local processing threads.

For online mode:

The useful options in online mode are:

   - <code>-t _num_</code> in online mode specifies the number of threads to use. It is up to the particular computation being run to allocate them as communication or local processing threads.
   - `-x` indicates that you wish to use an XOR-shared memory instead of the default additive-shared memory.  Some computations can work with either type (and respect this option) while others cannot (and ignore this option).

The args in online mode specify what computation to do.  These are currently mostly just unit tests.  The interesting one for now:

  - <code>duoram _depth_ _items_</code>: in a memory of size 2<sup>_depth_</sup>, do _items_ dependent updates, _items_ dependent reads, _items_ independent reads, _items_ independent updates, _items_ dependent writes, and _items_ dependent interleaved reads and writes.  Each batch is timed and measured independently.

It is vital that all three players are told the same computation to run!

The output of online mode (for each player) includes timings (real and CPU time), the number of messages sent, the number of message bytes sent in total, the Lamport clock (the number of network latencies experienced), and the number of local AES operations performed (each AES operation is an encryption of a single block), as well as the list of the number of precomputed resources used in each communication thread.

Example: In three terminals on the same host, run:

  - `./prac -p 0`
  - `./prac -p 1 localhost`
  - `./prac -p 2 localhost localhost r20:90 p:8`

to complete the preprocessing required for the following computation,
using 8 threads, and then:

  - `./prac -t 8 0 duoram 20 10`
  - `./prac -t 8 1 localhost duoram 20 10`
  - `./prac -t 8 2 localhost localhost duoram 20 10`

to do the online portion of the computation, also using 8 threads (which
will be local processing threads for the duoram computation).

For online-only mode:

The options and arguments are the same as for online mode, only you also add the `-o` option, and you don't have to have run preprocessing mode first.  At the end of the run, PRAC will output the number of resources of each type created, in a format suitable for passing as the arguments to preprocessing mode.

Example: In three terminals on the same host, run:

  - `./prac -o -t 8 0 duoram 20 10`
  - `./prac -o -t 8 1 localhost duoram 20 10`
  - `./prac -o -t 8 2 localhost localhost duoram 20 10`

No separate preprocessing step is needed.

## Docker instructions:

  - `cd docker`
  - `./build-docker`
  - `./start-docker`
    - This will start three dockers, each running one of the parties

Then to simulate network latency and capacity (optional):

  - `./set-networking 30ms 100mbit`

To turn that off again:

  - `./unset-networking`

If you have a NUMA machine, you might want to pin each player to one
NUMA node.  To do that, set these environment variables before running
`./run-experiment` below:

  - `export PRAC_NUMA_P0="numactl -N 1 -m 1"`
  - `export PRAC_NUMA_P1="numactl -N 2 -m 2"`
  - `export PRAC_NUMA_P2="numactl -N 3 -m 3"`

Adjust the numactl arguments to taste, of course, depending on your
machine's configuration.  Alternately, you can use things like `-C 0-7`
instead of `-N 1` to pin to specific cores, even on a non-NUMA machine.

Run experiments:

  - <code>./run-experiment _opts_ _args_</code>
    - opts and args are those of `./prac` above, *not including* the
      player number and other players' addresses, which are filled in
      automatically.

When you're all done:

  - `./stop-docker`
    - Remember this will lose anything you've saved on the docker filesystem, most notably any precomputed resources.
