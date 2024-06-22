# Synthesizing BC over unbounded domains via homogenization.

This file describes the research artifact for the publication:

On Completeness of SDP-Based Barrier Certificate Synthesis over Unbounded Domains
Hao Wu, Shenghua Feng, Ting Gan, Jie Wang, Bican Xia and Naijun Zhan
FM 2024 (Embedded System track)
Paper 113

## Preparation

- Mosek solver (v10.2) for solving SDPs.
- Julia (v1.9.2) for formulating Sum-of-Squares relaxations and translating them into SDPs
- Mathematica (v12) or Wolfram Engine (v14.0) for posterior verification and plotting graphs.

(Skip this part if you are using the virtual machine image.) Follow the instructions to prepare the environment:
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
5. Optional, if you want to verify the results and plot the graphs. Install Mathematica (commercial, https://www.wolfram.com/mathematica/) or Wolfram Engine (free but needs a developer licence, https://www.wolfram.com/engine/). 

## Virtual Machine

This artifact is distributed as a virtual machine (VM) in OVA format. The VM contains an installation of Ubuntu 22.04, along with the required software. The VM was created and tested using Oracle VirtualBox 7.0.18 running on a Windows 10 laptop with 16 GB of RAM and an intel i9-12900K CPU 2.50 GHz. During production the VM was allocated 4 GB of RAM and 1 processor through VirtualBox. The VM account is `vboxuser` and the password is `qwer12321`.

Open terminal and run `cd BCunbound`

## Important Files

```
BCunbound/
- Benchmarks/      % all benchmarks written in Julia
- Results/         % SDP computation results
   - system/       % system information for Mathematica to read
   - sound/        % results using Thm 3 
   - complete/     % results using Thm 5
   - completesemi/ % results using Thm 7
   - plots/        % trajectory plots using Mathematica
- run.jl           % solver SOS programs (Thm 3,5,7) using Mosek 
- verify.wls       % Molfram Engine script for verification
- verify.nb        % Mathematica notebook file for verifcation
- clean.sh         % script for clean results
```

## Evaluation Instructions
The purpose of the artifact is to substantiate the experimental results in Table 1 and Fig 1. The experiment consists of two parts: (1) Compute barrier certificate candidates within the degree range using Julia; (2) Verify the candidates using Mathematica/Wolfram Engine.   


- To compute barrier certificate candidates for all benchmarks, run `julia run.jl`. It will search for polynomial barrier certificates from degree 1 to 6 (using **Thm. 3** and **Thm. 5**) and semialgebraic barrier certificates from degree 1 to 4 (using **Thm. 7**). The returned candidates as well as the time overhead will be stored in `Results/` directory.
- To verify results in the first step, run `wolframscript verify.wls > result.txt`. The wolframe engine reads the results in `Results/` directory and tries to prove the barrier certificate candidates satisfy the barrier certificate conditions. If not, a point violating the conditions will be found by using `FindInstance` function. 

The output file `result.txt` will contain all the information within Table 1. Finally, the plots of all benchmarks can be found within the `Results/plots` folder.

## Claims
1. In our paper, all experiments were performed on a Mac lap-top with Apple M2 chip and 8GB memory. When preparing the artifact, we find that the verification step in VM typically takes more time. Though the verification time was not reported in our paper, this can lead to more question marks in Table 1 for benchmark `lorenz` and `lotka`. 

2. As explained in our paper, the time in Table 1 does not contain the verification time. Moreover, if a valid barrier certificate is found at degree 4, we do not count the SDP time overhead for degree 5 and 6.

3. Due to some unknown reasons, the Wolfram engine may crush after running for 20~30 minutes in VM (especially when verifying `lorenz` and `lotka` results). If this happens, please run `wolframscript verify_lorenz.wls > result_lorenz.txt` or `wolframscript verify_lotka.wls > result_lotka.txt` to separately verify the `lorenz` and `lotka` benchmarks.