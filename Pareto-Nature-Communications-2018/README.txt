Dynamic Programming for Computing the Pareto Frontier (for Energy and GHGs) on Tree Structured Networks  
Modified 6/17/2019
Version 1.0

================================================
Copyright 2019 Institute for Computational Sustainability

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

Jonathan M.Gomes-Selman(*) (Stanford University Dept. of Computer Science) jgs8@stanford.edu
Qinru Shi (Cornell University Dept. of Computer Science) qs63@cornell.edu
Guillaume Perez(*)  (Cornell University Dept. of Computer Science) guillaume.perez06@gmail.com 
(*) work performed while at Cornell

*** Citation ***

Please cite the following publications:

Jonathan M. Gomes-Selman, Qinru Shi, Yexiang Xue, Roosevelt Garcia-Villacorta, Alexander S. Flecker, Carla P. Gomes. Boosting Efficiency for Computing the Pareto Frontier on Tree Structured Networks. Proc. of conference on Constraint Programming, Artificial Intelligence and Operations Reserach 2018: 263-279

Xiaojian Wu, Jonathan Gomes-Selman, Qinru Shi, Yexiang Xue, Roosevelt Garcia-Villacorta, Elizabeth Anderson, Suresh Sethi, Scott Steinschneider, Alexander Flecker, Carla P. Gomes. Efficiently Approximating the Pareto Frontier: Hydropower Dam Placement in the Amazon Basin. AAAI 2018: 849-859 

Rafael M. Almeida, Qinru Shi, Jonathan M. Gomes-Selman, Xiaojian Wu, Yexiang
Xue, Hector Angarita, Nathan Barros, Bruce R. Forsberg, Roosevelt
Garcia-Villacorta, Stephen K. Hamilton, John M. Melack, Mariana
Montoya, Guillaume Perez, Suresh A. Sethi, Carla P. Gomes,
Alexander S. Flecker. Reducing greenhouse
gas emissions of Amazon hydropower with optimal dam planning, working
paper, 2019.

Please cite related works as appropriate in follow-up publications.

*** Compilation ***

Open a terminal window: The following instructions assume that the user will be issuing commands from a terminal window. 
 
There are two options for compiling and running the code
1) If g++ is installed correctly, run the following command within the source code directory:

   g++ -std=gnu++11 -pthread DP_Algorithm.cpp HyperNet.cpp main.cpp Pareto_Solution.cpp -o Amazon-E-GHG

2) Otherwise, make sure cmake and some c++ build toolchain/compiler (with c++11 support) are installed.
   Then run the two commands:

   cmake ./
   cmake --build ./

Both methods will generate an executable called "Amazon-E-GHG"

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

Move to the directory where the executable Amazon-E-GHG is. Then issue the following command:

./Amazon-E-GHG [-h] [OPTIONS] -criteria MAX/MIN MAX/MIN -basin BASIN_NAME path BASIN_FILENAME -epsilon ROUNDING_EPSILON

Options:
    Important notation
    ------------------------
    max/min = a string that is either "max" or "min"
    ------------------------

    -criteria <max/min> <max/min>
     (required) Reads in whether to maximize or minimize each criteria where the 
     order of the criteria matches that of the river basin file used as input.
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

./Amazon-E-GHG -criteria max min -batch 1000000 -basin Amazon -path AB_energy_ghg20.txt -epsilon 0.05 -thread 6

[the output for this example solution (.sol file) is provided in the downloaded zip folder. It can be easily read on a terminal window or a text editor] 
[epsilon = 0.1 (90% accuracy) should take around 13 seconds; epsilon = 0.05 (95% accuracy) should take around 22 seconds; epsilon = 0.01 (99% accuracy) should take around 97 seconds. If epsilon = 0 (exact version), the run time will be substantially longer (about 15 minutes) and the output file is much larger. These run times were obtained using a MacBook Pro 3.1 GHz Intel Core i5 16GB RAM]

NOTE: The order of the criteria matches that of the input file "AB_energy_ghg20.txt". For example, if in the the input
basin file the criteria are in the order "energy", "GHG" then the above command will run the algorithm to 
maximize energy and minimize greenhouse gas emisions(GHG).

*** River Basin File Format ***

Lines starting with % are comments
All fields are separated with a single space;

The file must contain 

1 - A line starting starting with "p" that specifies the problem in terms of:  number of nodes (where nodes represent connected river components), number of edges (where an edge represents a dam)

p #nodes #edges 

p 7 6


2 - A line starting with "dam_criteria" containing the name of the 2 criteria and status at the end (this indicates that the last column of each damn represents whether it is built or not built)

dam_criteria energy ghg status 


3 -  #dams lines starting with a "d" that specifies the information for each dam with the criteria in the order described in step (2): dam_id (number); energy; Greenhouse Gas emission; state (1:existing dam; 0: not existing yet). (Note: #dams has to be equal to the number of edges)
Example: 

d 1 30.0 10.0 1


4 - one line starting with "r" specifying the node# that corresponds to the root of the tree

r 0


5 - #edges lines starting with an "e" specifying, for each edge the parent-child nodes and the dam_id associated with the edge.(Note: #dams has to be equal to the number of edges)

e  parent-node# child-node# dam_id

e 1 2 1


Example:

File name Example.txt

p 7 6 2
dam_criteria energy ghg status
d 0 300 10 1
d 1  20 1 0
d 2 1000 5 0
d 3  10 1 0
d 4  100 1 0
d 5 100 1 0  
r 0
e 0 1 0
e 0 2 1
e 1 3 2
e 1 4 3
e 2 5 4
e 2 6 5

*** Output File Format ***
1) The date and time when the experiment was run

Date/time: Tue_Jun_18_15_45_17_2019


2) The river basin file used for the experiment

Data file: AB_energy_ghg20.txt


3) The wall time for the experiment (This is the true time of the experiment)

Wall time: 3.6614 seconds.


4) The cpu time of the experiment. This captures how much "time" the cpu used, 
   which factors in time on multiple threads.

CPU time: 11.9694 seconds.


5) The random seed used for the experiment

seed: 1560887114


6) The number of Pareto Optimal solutions found

num_solutions: 55


7) The number of nodes in the tree, which represents the number
   of times that pruning occurred

# pruning steps (# nodes): 608


8) How many possible partial solutions are possible in the tree
   considering all the dam combinations

Max policies considered: 1508652


9) The number of partial solutions actually considered. Note when
   using a non-zero epsilon (rounding based approximation) this number
   can be smaller than the number of Max policies considered because
   rounding causes solutions to collide/become equal.

Policies considered: 3883616


10) The number of partial solutions that were pruned by the algorithm

Pruned policies: 3758926


11) The epsilon used for the approximation. Note, if epsilon = 0 then
    we compute the exact Pareto Frontier.

epsilon: 0.25

12) The batch size used for pruning in the Dynamic Programming algorithm

batch size: 100000


13) The criteria considered in the experiment, listed in the order that
    they appear in the generated solutions

criteria: energy, ghg

14) The remaining N lines (where N is equal to the number of solutions generated,
    given on line 6) include the solutions found. Each solution is described in 
    comma separated format as follows: 

(first criteria value), (second criteria value), (number of dams built), [list of dam ids built]

Example:
42876.6, 1190.65, 268, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 130 131 139 143 151 152 156 157 158 161 162 165 167 185 186 188 189 190 196 200 207 222 227 230 232 233 236 237 238 241 242 246 249 251 252 253 254 255 257 262 269 270 271 274 284 291 292 293 294 295 296 302 304 306 309 310 312 313 315 318 320 321 323 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 352 353 354 355 356 357 358 359 361 362 364 365 366 370 371 372 373 376 377 378 379 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 421 431 435 444 453 468 470 472 476 477 478 486 487 488 489 491 493 497 498 499 500 501 504 505 506 507 508 509

*** Main contacts ***

If you have any questions or suggestions, please feel free to contact JGS at jgs8@stanford.edu or QS at qs63@cornell.edu or gomes@cs.cornell.edu
