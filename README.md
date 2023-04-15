# ColorfulStarLocal

> This repository is the official implementation of ***Parallel Colorful ℎ-star Core Maintenance in Dynamic Graphs***

> The default graph coloring algorithm: **Degree**, i.e. following a non-increasing order of node degrees

## Getting the code
You can download a copy of all the files in this repository by cloning the [git](https://git-scm.com/) repository:
```
https://github.com/SGAOCR/ColorfulStarLocal.git
```
or [download a zip archive](https://github.com/Gawssin/ColorfulStarLocal/archive/refs/heads/main.zip).

## Dependencies
```
gcc: 11.2.0 or above
OpenMP: 4.5 or above
```

## Repo Structure
- **Algorithms**
  - *The implementation of different algorithms*
- **header**
  - *Header Files*

## Algorithm Description
We here consider 7 different algorithms for the colorful h-star core decomposition and maintenance.
| Algorithm  |  Shorthand    | Description | 
| ----------- | ----------- | ----------- |
| The Peeling method     |  **corePeel** | The basic method to compute the colorful h-star core decomposition following the peeling paradigm, that is to remove the node with minimum colorful h-star degree   |
| The local algorithm | **Local**    | Compute the colorful h-star core decomposition by iteratively computing the colorful h-star H-index for each node until all of them converge.       |
| Optimization: Asynchronous computing  | **OPT-1**    | To compute the colorful h-star H-index of node u, allways use the latest H-indexes of its neighbors.        |
| Optimization: Processing ordering heuristic  | **OPT-2**    | In each iteration,  process nodes following an increasing order of their (𝑖 − 1)-order H-indexes instead of the random processing order.        |
| Optimization: Pruning technique  | **OPT-3**    | Detect redundant computation and skip nodes whose 𝑛-order H-indexes will not be changed in each iteration. |
| The edge deletion algorithm  | **EdgeDel**    | efficiently update colorful h-star core numbers of affected nodes after a deletion of an edge. |
| The edge insertion algorithm  | **EdgeIns**    | efficiently update colorful h-star core numbers of affected nodes after an insertion of an edge. |


## Reproducibility

### 1. The Peeling method
#### To compile
```
$ cd ./Algorithms/
$ g++ -O3 -fopenmp -o corePeel corePeel.cpp
```
#### To Run
```
$ ./corePeel h FilePath colorAlgo [basic]
```
| Argument  | Description |
| :-----| :---- |
| **corePeel** | executable file |
| **h** | the size of stars |
| **FilePath** | input file path |
| **colorAlgo** | graph coloring algorithm: <br>"0" for **Degree**; <br>"1" for **Degen**; <br>"2" for **FF**; <br>"3" for **SD**. <br> Find more details [here](https://ieeexplore.ieee.org/document/9835239/). |
| **[basic]** | *Optional*, indicate which algorithm will be performed. <br>"basic" for the basic version **HStarDP**, i.e. recomputing the colorful h-star degree when a neighbor is removed; <br>the default is the advanced version **HStarCD**, i.e. using the proposed updating technique to calculate the colorful h-star degree. |

### 2. The local algorithm
#### To compile
```
$ cd ./Algorithms/
$ g++ -O3 -fopenmp -o Local Local.cpp
```
#### To Run
```
$ ./Local h FilePath N_threads
```
| Argument  | Description |
| :-----| :---- |
| **Local** | executable file |
| **h** | the size of stars |
| **FilePath** | input file path |
| **N_threads** | the number of threads |

### 3. Three optimizatons
#### To compile
```
$ cd ./Algorithms/
$ g++ -O3 -fopenmp -o LocalOPT LocalOPT.cpp
```
#### To Run
```
$ ./LocalOPT h FilePath N_threads OPTtype
```
| Argument  | Description |
| :-----| :---- |
| **LocalOPT** | executable file |
| **h** | the size of stars |
| **FilePath** | input file path |
| **N_threads** | the number of threads |
| **OPTtype** | indicate which optimization will be applied: <br>"1" for **OPT-1**; <br>"2" for **OPT-1 + OPT-2**; <br>"3" for **OPT-1 + OPT-2 + OPT-3** (also stands for OPT*).   |


### 4. The edge deletion algorithm
#### To compile
```
$ cd ./Algorithms/
$ g++ -O3 -fopenmp -o EdgeDel EdgeDel.cpp
```
#### To Run
```
$ ./EdgeDel h FilePath N_Del RecompFlag N_threads
```
| Argument  | Description |
| :-----| :---- |
| **EdgeDel** | executable file |
| **h** | the size of stars |
| **FilePath** | input file path |
| **N_Del** | the number of edges to be deleted |
| **RecompFlag** | indicate whether to  compare with ("1") the recomputing method, i.e., call **Local** compute for all nodes or not ("0")|
| **N_threads** | the number of threads |


### 5. The edge insertion algorithm
#### To compile
```
$ cd ./Algorithms/
$ g++ -O3 -fopenmp -o EdgeIns EdgeIns.cpp
```
#### To Run
```
$ ./EdgeIns h FilePath N_Ins RecompFlag N_threads
```
| Argument  | Description |
| :-----| :---- |
| **EdgeIns** | executable file |
| **h** | the size of stars |
| **FilePath** | input file path |
| **N_Ins** | the number of edges to be inserted |
| **RecompFlag** | indicate whether to  compare with ("1") the recomputing method, i.e., call **Local** compute for all nodes or not ("0")|
| **N_threads** | the number of threads |
