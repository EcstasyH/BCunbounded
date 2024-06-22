# On Completeness of SDP-Based Barrier Certificate Synthesis over Unbounded Domains
Synthesizing BC over unbounded domains via homogenization. 

## Preparation 
- Mosek solver (v10.2) for solving SDPs.
- Julia (v1.9.2) for formulating Sum-of-Squares relaxations and translating them into SDPs
- Mathematica (v12) or Wolfram Engine (v14.0) for posterior verification and plotting graphs.

Either use the virtual machine or follow the instructions to prepare the environment:
1. Install Mosek solver (https://www.mosek.com/downloads/) and apply for an academic license via (https://www.mosek.com/products/academic-licenses/)
2. Install Julia programming language (https://julialang.org/)
3. Clone this repository from github
   - run `git clone https://github.com/EcstasyH/BCunbounded`
   - run `cd BCunbounded` 
4. Install necessary julia packages
   - run `julia` to start an interactive session 
   - run `import Pkg`
   - run `Pkg.activate(".")` to activate the environment
   - run `Pkg.instantiate()` to download necessary dependencies
   - run `exit()` to quit the interactive session
5. (optional, if you want to verify the results and plot the graphs) install Mathematica
   - Mathematica becomes inefficient when the problem instance involves too many variables and/or the degree is too large. 

## Experiments
The experiments consist of two parts: (1) Compute barrier certificate candidates within the degree range using Julia (2) verify the candidates using Mathematica.   

- To compute barrier certificate candidates for all benchmarks, run `julia run.jl`.
- To verify results in the first step, run `wolframscript verify.wls > result.txt`


## Verify Results
The results for benchmarks all stored in `Results` directory. We use Mathematica to read these results and verify them. Please refer to the content of `Verify.nb` file. 

## File Structure