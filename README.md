# Porechop_ABI

Porechop_ABI (*ab initio*) is an extension of [Porechop](README_PORECHOP.md) whose purpose is to process adapter sequences in ONT reads. 

The difference with the initial version of Porechop is that Porechop_ABI does not use any external knowledge or database for the adapters.  Adapters are discovered directly from the reads using approximate k-mers counting and assembly. Then these sequences can be used for trimming, using all standard Porechop options.

The software is able to report a combination of distinct sequences if a mix of adapters is used. It can also be used to check whether a dataset has already been trimmed out or not, or to find leftover adapters in datasets that have been previously processed with Guppy.

Note that Porechop_ABI is not designed to handle barcoded sequences adapters.
Demultiplexing should be done using standard Porechop commands or other appropriate tools.  

# Table of Content
* [Requirements](#requirements)
   * [Operating System version](#operating-system-version)
   * [Rust (only to build from source)](#rust-only-to-build-from-source)
   * [Python version](#python-version)
   * [Networkx](#networkx)
* [Installation](#installation)
* [Quick usage](#quick-usage)
* [Advanced usage](#advanced-usage)
   * [General purpose options](#general-purpose-options)
   * [Algorithm options](#algorithm-options)
* [Config File](#config-file)
   * [List of possible parameters](#list-of-possible-parameters)
* [Contributors](#contributors)
* [References](#references)
* [License](#license)


## Requirements 

If you want to install Porechop_ABI, your system must satisfy some requirements. Porechop_ABI is written in Python and Rust. The Rust part (the `porechop_abi._native` extension module, sources in `rust/`) does the sequence alignment, the approximate k-mer counting and the multiple alignment consensus; networkX is used for graph assembly.
In the case you chose to install using Conda, you can ignore this section (except the [OS](operating_system_version) requirement) and just use an up-to-date version of conda to install Porechop_ABI.
Check the [installation](INSTALL.md) guide for more details.


### Operating System version
 
* Linux, any version with compatible python requirements should work. 
* Mac OS >= 10.12.


### Rust (only to build from source)
Building from source compiles the extension module, which requires a Rust toolchain, version 1.85 or later (install it with [rustup](https://rustup.rs), check with `cargo --version`).
No C++ compiler, SeqAn or zlib is needed. A built wheel contains the compiled module and needs no Rust to install.

It reproduces the results of the former C++ components (SeqAn 2.4 alignment, `approx_counter` and `msa_consensus`) exactly, and is much faster.


### Python Version
Porechop is written mainly in Python3.
Python version must be 3.8 or above.<br>

With [X] being the version of python you want to install:

On Linux
~~~
sudo apt-get install python3.[X]
~~~

On Mac os (under macports):
~~~
sudo port install python3.[X]
~~~

### Networkx

With Pip:
~~~
pip install networkx
~~~
If you have permission errors, you may want to install it for the current user only:
~~~
pip install --user networkx 
~~~

On linux:
~~~
sudo apt-get install python-networkx
~~~
On Mac os (under macports):
~~~
sudo port install py3[X]-networkx 
~~~
With [X] being the version of python you use.


## Installation 

Porechop_ABI can be installed by several means. 
The easiest way is to use the conda package management system and install it from the bioconda channel in a clean environment.
All details about installation procedures are listed in the [installation](INSTALL.md) guide.<br>


## Quick usage 


### Adapter inference and trimming

`porechop_abi --ab_initio -i input_reads.fastq -o output_reads.fastq` <br>
`porechop_abi -abi -i input_reads.fastq -o output_reads.fastq`


This command-line with the -abi flag allows to first guess the adapters from the reads, add the adapters to the list of Porechop adapters  (adapters.py file) and then run Porechop as usual. 
 It is compatible with all Porechop options for trimming. The input_reads.fastq file shoud contain the set of raw ONT reads. The resulting trimmed reads are saved in  output_reads.fastq file.


### Adapter inference only 

`porechop_abi -abi --guess_adapter_only -i input_reads.fastq.gz -o output_reads.fastq`
`porechop_abi -abi -go -i input_reads.fastq -v 0 -o output_reads.fastq.gz`

This command-line, with the additional -go flag,  allows to only guess and print the adapter sequences from the reads. It then stops the execution of the program, without 
trimming the reads. 

### Test files

Two read  files are provided to test Porechop_ABI.


__Faster test:__ Simulated data (with reduced sampling, discarding default database) 

`porechop_abi -abi -go -ddb -nr 5 -cr 15 -i test/test_simulated_10k_read.fasta -tmp /tmp/pabi_temp -o /dev/null`

__Slower test:__ Real data (with standard parameters, discarding default database) 

`porechop_abi -abi -go -ddb -i test/test_realdata_10k_read.fasta -tmp /tmp/pabi_temp -o /dev/null`

[//]: # (TODO: Add expected results for P_ABI tests.)


### Warnings


Porechop_ABI can output two types of warning during the execution of adapter inference.

Low frequency warning: This warning means that the k-mers used to build the adapter sequences do not show a clear over-representation. In this case, the result of the algorithm is questionable, because the signal is not reliable. This may happen when the input reads are already trimmed, for example. 

Poor consensus warning: This warning is activated when Porechop_ABI finds more than two distinct start (or end) adapters and when each of them is found in less than 30% reads. 



## Advanced usage
Porechop_ABI allows to tune parameters, either for the ABI phase (ab initio) or the read trimming phase. We decribe hereafter the options that are specific to the ABI phase. Some of them are accessible in the config file, and some of them in the command-line.
For all usages and descriptions regarding the read trimming phase, you can refer to the Porechop [documentation](README_PORECHOP.md).

### General-purpose options 

#### Path to the config file  
  
```bash
-abc / --ab_initio_config [path_to_file]
```
Allows you to set a custom config file for the ab initio phase (default file in Porechop folder).
The config file comes with its own set of parameters. Check the "Config File" section for more details.

#### Temporary files

```bash
-tmp / --temp_dir [path_to_folder]
```
Path to a writable temporary directory, used to store count file and temporary fasta files.
The provided directory will be created if it does not exists, if the path is writable.
By default, Porechop_ABI creates a ./tmp/ folder in the working directory (where the command is launched).

#### Random seed

```bash
--seed [integer]
```
Seed for the random sampling of reads in the ab initio phase. Two runs with the same seed and the same input give the same result.
By default, the sampling is random.

#### Custom adapters

```bash
--custom_adapters / -cap [path_to_file]
```     

By default, Porechop uses a database of adapters, stored in the adapter.py file. The flag -cap  allows you to use a  set custom adapters without editing the adapter.py file.
These adapters are only used at runtime, not added to the adapter.py file.
Custom adapters must be stored in a text file following this format:
```    
    line 1: Adapter name
    line 2: Start adapter sequence
    line 3: End adapter sequence
    --- repeat for each adapter pair---
```

If your adapters do not contain the start or end sequence,
just put an empty line.

Each adapter will be named using the adapter name as a prefix, like:

```
"adapter_name_Top"
"adapter_name_Bottom"
```

#### discarding Porechop native databse

```bash
--discard_database  / -ddb
```                              
In case you supply your own adapters, it can be useful (and faster) to ignore the adapters from the Porechop database.
This option was added for this situation, and require either ab initio (-abi) or a custom adapter (-cap) to be set. Default: False


### Algorithm options 

The algorithm implemented in the ABI module depends on a series of parameters concerning k-mer selection, sampling, and assembly. They all have default values that work well in practice.

* Number of reads in each sample. By default, the value is 40,000 reads.
* Length of start and end regions selected for each read. By default, the value is 100 nucleotides.
* Length of the k-mers. By default, the value is 16.
* Number of top k-mers in the selection of frequent k-mers. By default, the value is 500.
* Low complexity threshold: The value depends on the k-mer length. By default, it equals 1 for k-mers of length 16. It is adjusted automatically for other values of k.
    
All these values can be customized by the user in the [Config file](#config-file). 
Additional parameters can be specified via the command-line.

```bash
--export_graph [output_file]
```    
Let you export the assembly graph with some annotations to the provided path.
Graphs are exported in graphml format, and can be visualised with programs such as Gephi.


*The following options are specific to the drop-cut algorithm used to adjust sequences length after the (greedy) assembly step.
Only experienced users shoud try to adjust these settings.*

```bash
--window_size / -ws  [int]
```   
Change the size of the smoothing window used in the drop cut algorithm. (set to 1 to disable smoothing).


```bash
--no_drop_cut / -ndc
```   
Disable the drop cut step entirely (default: False)


```bash
--number_of_run / -nr  [int]
```   
Porechop_ABI builds consensus sequences from multiple samples to increase confidence in the inferred adapter sequence. The -nr flag allows you to specify the number of initial samples. Each count file is exported separately, and can be reviewed if needed. Default: 10, set to 1 for single sample mode.

```bash
--consensus_run / -cr [int]
```
If using multi-run option above 1, set the number of additional samples performed if no stable consensus is  found after the first trial. Default: 20.

After this step, all inferred adapter sequences are clustered by similarity, and consensus are generated for each cluster.


```bash
--export_consensus / -ec [output_file]
```
Path to export the intermediate adapters found in consensus mode.


```bash
--all_above_x / -aax [int]
```
Only select consensus sequences if they are built using at least x percent of the total adapters.
Default : 10 (3 sequences out of  30). You can set it to 0 to keep all consensus.
  
```bash
--best_of_x / -box [int]
```
Only select the best x consensus sequences from all consensus found. Default: 0, reports all consensus.


## Config file

As seen before, some advanced options 
are accessible in the config file located at /Porechop_ABI/porechop/ab_initio.config. Note that default values work just fine in practice, and it is most likely that you will not need to edit this file.

If you installed Porechop_ABI with pip, take note that the actual config file used when using porechop from command line will be the one stored in your installation folder (example: /usr/local/lib/python3.__X__/dist-packages/porechop/porechop/ab_initio.config).

In the case you still want to use custom options, it is recommanded to pass a config file as a parameter to Porechop_ABI using -abc / --ab_initio_config instead of editing the default config.


```bash
porechop_abi --ab_initio_config YOUR_CONFIG.config ...
```
### List of possible parameters:
```
  // Sampling options
  sn        : Number of read to sample (default: 40 000 reads)
  sl        : Size of the samples (default : 100 bases)
  limit     : Number of k-mer kept after counting. (default: 500)
  solid_km  : Use solid kmers instead of a hard limit (specify the minimum count for a k-mer, will override "limit" option).
  
  // K-mer options
  k           : K-mer size (between 2 and 32, default 16).
  lc          : Low Complexity filter value (default 1.0).
  forbid_kmer : Forbidden k-mer list. (list of k-mers, one per line).
  exact_out   : Also export the exact count (mainly for debug / control)
  
  
  // Misc
  v           : Verbosity level (0: minimal output, 1: normal, 2: debug).
  nb_thread   : Number of threads used for the approximate counting.
  skip_end    : Do not process end adapters. (mostly for debug)
  nb_of_runs  : Number of time the sampling and count is performed. (Default: 10)
```

**Notes:**
Verbosity, number of threads, and number of runs are copied from the arguments supplied to Porechop_ABI.
Changing the value of these parameters in the config file will have no effect.
If you need to adjust one of them, refer to porechop -h


## Contributors

The original Porechop program was provided by Ryan Wick.
The ab initio extension is developed by Quentin Bonenfant, Laurent Noé and Hélène Touzet.


## References
If you want more details about the algorithm, internal working and performances of this tool, a preprint is available here:<br>
https://www.biorxiv.org/content/10.1101/2022.07.07.499093

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
