Computing the exact and approximate Pareto frontier on tree-structured networks with application to reducing the adverse impacts of hydropower expansion on ecosystem services in the Amazon Basin  
Modified 11/1/2021
Version 1.1

================================================
Copyright 2021 Institute for Computational Sustainability

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
================================================

*** Main Contributors ***

Jonathan M.Gomes-Selman(*) (Stanford University Dept. of Computer Science) jgs8@cs.stanford.edu
Qinru Shi (Cornell University Dept. of Computer Science) qs63@cornell.edu
Guillaume Perez(*)  (Cornell University Dept. of Computer Science) guillaume.perez06@gmail.com 
(*) work performed while at Cornell

*** Abstract***

Multi-objective optimization plays a key role in the study of real-world problems, as they often involve multiple criteria. In multi-objective optimization, it is important to identify the so-called Pareto frontier, which characterizes the trade-offs between the objectives of different solutions. We provide a C++ implementation of exact and approximate dynamic programming (DP) algorithms for computing the Pareto frontier on tree-structured networks. The code uses a specialized divide-and-conquer approach for the pruning of dominated solutions. This optimization outperforms the previous approaches, leading to speed-ups of two to three orders of magnitude in practice.

We apply a rounding technique to the exact dynamic programming algorithm that provides a fully polynomial-time approximation scheme (FPTAS). The FPTAS finds a solution set of polynomial-size, which approximates the Pareto frontier within an arbitrary small e factor and runs in time that is polynomial in the size of the instance and 1/ep.  

We illustrate the code by evaluating trade-offs in ecosystem services due to the proliferation of hydropower dams throughout the Amazon basin. In particular, we apply our algorithms to identify portfolios of hydropower dam sites that simultaneously minimize impacts on river flow, river connectivity, sediment transport, fish biodiversity, and greenhouse gas emissions while achieving energy production goals, at different scales, including the entire Amazon basin. The code can be easily adapted to compute the Pareto frontier of various multi-objective problems for other river basins or other tree-structured networks. This work is described in the manuscript  by Flecker et al., entitled “Reducing adverse impacts of Amazon hydropower expansion” in press, Science, 2021.

*** Citation ***

Please cite the following publications:

Jonathan M. Gomes-Selman, Qinru Shi, Yexiang Xue, Roosevelt Garcia-Villacorta, Alexander S. Flecker, Carla P. Gomes. Boosting Efficiency for Computing the Pareto Frontier on Tree Structured Networks. Proc. of conference on Constraint Programming, Artificial Intelligence and Operations Reserach 2018: 263-279

Xiaojian Wu, Jonathan Gomes-Selman, Qinru Shi, Yexiang Xue, Roosevelt Garcia-Villacorta, Elizabeth Anderson, Suresh Sethi, Scott Steinschneider, Alexander Flecker, Carla P. Gomes. Efficiently Approximating the Pareto Frontier: Hydropower Dam Placement in the Amazon Basin. AAAI 2018: 849-859 

Alexander S. Flecker, Qinru Shi, Rafael M. Almeida, Héctor Angarita, Jonathan M. Gomes-Selman, Roosevelt García-Villacorta, Suresh A. Sethi, Steven A. Thomas, N. LeRoy Poff, Bruce R. Forsberg, Sebastian A. Heilpern, Stephen K. Hamilton, Jorge D. Abad, Elizabeth P. Anderson, Nathan Barros, Isabel Carolina Bernal, Richard Bernstein, Carlos M. Cañas, Olivier Dangles, Andrea C. Encalada, Ayan S. Fleischmann,  Michael Goulding, Jonathan Higgins, Céline Jezequel, Erin I. Larson, Peter B. McIntyre, John M. Melack, Mariana Montoya, Thierry Oberdorff, Rodrigo Paiva, Guillaume Perez, Brendan H. Rappazzo, Scott Steinschneider, Sandra Torres, Mariana Varese, M. Todd Walter, Xiaojian Wu, Yexiang Xue, Xavier E. Zapata-Ríos, and Carla P. Gomes. Reducing adverse impacts of Amazon hydropower expansion. Science, 375(6582):753–760, 2022. 

Please cite related works as appropriate in follow-up publications.

*** Compilation ***

Open a terminal window: The following instructions assume that the user will be issuing commands from a terminal window. 
 
There are two options for compiling and running the code
1) If g++ is installed correctly, run the following command within the source code directory:

   g++ -std=gnu++11 DP_Algorithm.cpp Frontier_List.cpp HyperNet.cpp main.cpp Pareto_Opt_Node.cpp Pareto_Opt_List.cpp Pareto_Solution.cpp -o Amazon

2) Otherwise, make sure cmake and some c++ build toolchain/compiler (with c++11 support) are installed.
   Then run the two commands:

   cmake ./
   cmake --build ./

Both methods will generate an executable called "Amazon"

Notes on installation:
- CMake can be downloaded at https://cmake.org/download/. After downloading CMake, make sure to 
  install CMake for command line use. On a Mac, this can be done by:
      a) Opening the CMake application
      b) Select Tools from the top menu
      c) Select "How to Install for Command Line Use"
      d) Follow the instructions for the second option "Or, to install symlinks to '/usr/local/bin', run:"
         where you must copy and run the command "sudo ... --install" into a terminal window.

- For a mac, one quick way to ensure that you have a proper c++ compiler is to install the commandline developer tools (or all of Xcode) by running: 

    xcode-select --install

*** Usage ***

Move to the directory where the executable Amazon is. Then issue the following command:

./Amazon -criteria NUMBER_OF_CRITERIA CRITERIA_NAME ... CRITERIA_NAME -batch BATCH_SIZE -basin BASIN_NAME path INPUT_FILENAME -epsilon ROUNDING_EPSILON

Options:
    Important notation
    ------------------------
    max/min = a string that is either "max" or "min"
    ------------------------

    -criteria <int> <string> ... <string>
     (required) Reads in relevant criteria. The number of strings provided must equal the 
     first integer argument.
    -path <string> 
     (required) The file path with the data for the river basin
    -basin <string>
     (required) The name of the river basin that we are working on
    -seed <int> 
     The random seed to use in the experiment (NOTE: if this is not given a random seed is chosen)
    -epsilon <double>
     Real valued epsilon in the range [0.0, 1.0] used for the epsilon 
     approximation (rounding) in the experiments. For two criteria problems, 
     use epsilon=0 to get the exact Pareto frontier. For faster results, consider 
     using epsilon=0.05 (i.e. within 5% theoretical guarantee of the true Pareto Frontier).    
    -batch <int>
     Number of solutions to keep in a batch when pruning dominated partial solution. 
     Commonly used batch sizes are 1000000 or 10000000 depending on available memory.
    -thread <int>
     Number of threads to use for parallelization. A good default value is 6 
    
*** Example of Running the Executable ***

        ./Amazon -criteria 2 energy connectivity -basin Amazon -path AB_cppinput_2021_Feb.txt -epsilon 0.001 -batch 10000000
        ./Amazon -criteria 2 energy connectivity -basin Amazon -path AB_cppinput_2021_Feb.txt -epsilon 0 -batch 10000000
        ./Amazon -criteria 3 energy connectivity sediment -basin Amazon -path AB_cppinput_2021_Feb.txt -epsilon 0.1 -batch 10000000
	./Amazon -criteria 6 energy connectivity sediment biodiversity dor ghg -basin Napo -path napo_cppinput_2021_Feb.txt -epsilon 0.05 -batch 10000000
	./Amazon -criteria 6 energy connectivity sediment biodiversity dor ghg -basin Napo -path napo_cppinput_2021_Feb.txt -epsilon 0 -batch 10000000

[Note that we always include energy as the first criterion. The reason is that energy is the only criterion we have that encourages building more dams. If we do not include energy as a criterion, we get a single optimal solution where we build no dams.]
[The output files can be easily read on a terminal window or a text editor] 
[Optimizing for energy and connectivity for epsilon = 0.001 for the whole Amazon (99.9% accuracy) should take around 10 seconds; optimizing for energy and connectivity for epsilon = 0 (exact) should take around 2 minutes; optimizing for three criteria: energy, connectivity, and sediment for epsilon = 0.1 (90% accuracy) should take around 80 seconds; optimizing for six criteria for the Napo basin (a basin with 24 dams) for epsilon = 0.05 (95% accuracy) should take below 1 second; optimizing for six criteria for the Napo basin (a basin with 24 dams) for epsilon = 0 (exact) should take around 2 seconds.  These run times were obtained using a MacBook Pro 2.3 GHz Intel Core i5 16GB RAM. 
[* Currently we can only run for six criteria within reasonable time and error bound on smaller basins such as the Maranon or the Napo. DO NOT attempt to run for six criteria on the whole Amazon. For more information, refer to the last section of the README.]

*** Input File Format ***

Lines starting with % are comments
All fields are separated with a single space;

The file must contain 

1 - A line starting with "p" that specifies the problem in terms of:  the number of nodes (where nodes represent connected river components), the number of edges (where an edge represents a dam), the number of node criteria, the number of dam criteria. 

p #nodes #edges #node_criteria #dam criteria

Example:
p 510 509 2 6


2 - A line starting with "dam_criteria" containing the name of the criteria related to dams(the status criteria indicates whether a dam is already built or not). The code currently only accept criteria that are already included in the provided input files.
Example:

dam_criteria energy sed_trap status biodiversity ghg dor 


3 - A line starting with "node_criteria" containing the name of the criteria related to nodes. The code currently only accept criteria that are already included in the provided input files.
Example:

node_criteria connectivity sediment


4 -  #dams lines starting with a "d" that specifies the information for each dam. The first number after "d" is the index of the dam (dam_id). The rest of the numbers are values of the dam criteria in the order described in step (2).
Example: 

d 40 740.0 0.872701517 0 2.3811145017 13.4 94.95341541


5 -  #nodes lines starting with an "n" that specifies the information for each node. The first number after "n" is the index of the node. The rest of the numbers are values of the node criteria in the order described in step (3).
Example: 

n 0 30128378.7642 55.336082413

6 - one line starting with "r" specifying the node# that corresponds to the root node of the tree
Example:

r 0

7 - #edges lines starting with an "e" specifying, for each edge, the parent-child nodes and the dam_id associated with the edge. (Note: #dams has to be equal to the number of edges)

e  parent-node# child-node# dam_id

Example:

e 0 1 107


Example input file:

File name: example_input.txt

p 7 6 2 4
dam_criteria energy sed_trap status ghg
node_criteria connectivity sediment
d 1 300.0 0.5 1 1.1 
d 2  20.0 0.1 0 0.2
d 3 1000.0 0.6 0 14.5
d 4  10.0 0.0 0 0.0
d 5  100.0 0.2 0 1.5
d 6 100.0 0.3 0 0.3 
n 0 143.5 23.4
n 1 5.6 0.5
n 2 74.0 9.8
n 3 41.5 6.2
n 4 9.6 0.0
n 5 30.0 4.3
n 6 95.0 20.6
r 0
e 0 1 1
e 0 2 2
e 1 3 3
e 1 4 4
e 2 5 5
e 2 6 6

*** Output File Format ***
1) The date and time when the experiment was run

Date/time: Sat_Mar_13_15_39_10_2021

2) The river basin file used for the experiment

Data file: AB_cppinput_2021_Feb.txt


3) The wall time for the experiment (This is the true time of the experiment)

Wall time: 79.3914 seconds.

4) The cpu time of the experiment, which captures how much "time" the cpu used, 
   which factors in time on multiple threads.

CPU time: 149.1856 seconds.

5) The random seed used for the experiment

seed: 1615667870

6) The number of Pareto Optimal solutions found

num_solutions: 3221


7) The number of nodes in the equivalent binary tree, which represents the number
   of times that pruning occurred. (Should be slightly larger than the number of dams)

# pruning steps (# nodes): 606

8) How many possible partial solutions are possible in the tree
   considering all the dam combinations

Max policies considered: 25690824

9) The number of partial solutions actually considered. Note when
   using a non-zero epsilon (rounding based approximation), this number
   can be smaller than the number of Max policies considered because
   rounding causes solutions to collide/become equal.

Policies considered: 31928400

10) The number of partial solutions that were pruned by the algorithm

Pruned policies: 31689695

11) The epsilon used for the approximation. Note, if epsilon = 0 then
    we compute the exact Pareto Frontier.

epsilon: 0.1

12) The batch size used for pruning in the Dynamic Programming algorithm

batch size: 10000000

13) The criteria considered in the experiment, listed in the order that
    they appear in the generated solutions

criteria: energy, connectivity, sediment

14) The remaining N lines (where N is equal to the number of solutions generated,
    given on line 6) include the solutions found. Each solution is described in 
    comma separated format as follows: 

(first criteria value), (second criteria value), ..., (last criteria value), (number of dams built), [list of dam ids built]

Example:
43365.9, 3.66885e+07, 564.044, 321, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 112 113 114 116 117 118 122 128 133 134 140 142 143 144 148 149 151 156 162 164 165 167 171 179 180 181 183 184 185 186 187 188 189 190 192 193 194 195 196 197 198 199 200 202 204 205 206 207 208 211 212 214 222 223 224 230 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 269 270 271 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 299 300 301 302 303 304 307 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 348 349 350 356 363 365 370 372 373 375 376 378 379 382 383 384 385 386 387 388 389 390 391 392 393 394 396 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 438 439 444 453 455 463 470 472 473 474 481 498 499 504 505 506 507 508 509



*** Computing the Pareto Frontier for six criteria***
The runtime of this algorithm we developed is polynomial with regard to the number of dams but still exponential with regard to the number of criteria, which means that both the runtime and the number of solutions will increase dramatically as the number of objectives goes up. For instance, for two criteria (energy and ghg), we are able to compute the exact Pareto frontier (epsilon=0) for the whole Amazon within 20 minutes (wall-clock time, 36 threads; around 10 hours CPU time) depending on the criteria we choose; for five criteria, we are able to run the algorithm for (epsilon=0.4) for the whole Amazon in 17 hours (wall-clock time, 36 threads; around 9.3 days CPU time); 

When optimizing for 6 criteria for the Amazon basin, however, we can only run for large error margins like epsilon=1.5 or epsilon=2.0, and the runtimes are 2 days and 7 hours (wall-clock time, 36 threads), respectively. The solution we employed to deal with the heavy costs and large error margins associated with optimizing for more criteria is to complement the approximate Pareto frontier with optimization results for all possible 2 criteria and 3 criteria combinations. Note that when we do not include energy as a criterion, we get a single optimal solution where we build no dams. Therefore, we only need to run the algorithm for criteria combinations that include energy, which are five 2 criteria combinations and ten 3 criteria combination.


*** Main contacts ***

If you have any questions or suggestions, please feel free to contact JGS at jgs8@cs.stanford.edu or QS at qs63@cornell.edu or gomes@cs.cornell.edu


*** Github link ***

Further updates of the code can be found in the github repository of this project: 
https://github.com/gomes-lab/Pareto-Frontier
