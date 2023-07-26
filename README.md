
DISCLAIMER: AleRax is still under development and not published! Some features are already implemented, but the code is not stable. Please contact me if you want to start to work with AleRax.

# AleRax  

AleRax is a parallel tool for species tree - gene tree reconciliation under gene duplication, loss, and HGT. For each gene family, it takes as input a gene family tree distribution (typically inferred with Bayesian inference tools such as MrBayes, PhyloBayes, etc.). AleRax can perform the following operations:
* Maximum likelihood species tree inference
* Species tree rooting (using the DTL events)
* Model parameter (e.g. DTL event probabilities) inference
* Reconciled gene tree sampling (we reimplemented [`ALE`](https://github.com/ssolo/ALE)) 

We are also working on the following features:
* Relative order of speciation event (relative dating) from HGT constraints
* Highway of transfer inference

## Requirement

* A Linux or MacOS environnement
* gcc 5.0 or > 
* CMake 3.6 or >
* MPI (optional, but required if you want to use parallelization)

## Installation 


To download AleRax, please use git,  and clone with --recursive!!!

```
git clone --recursive https://github.com/BenoitMorel/AleRax
```

For some reason, cloning coraxlib (one of the dependencies) sometimes fails:
```
Cloning into '/home/benoit/github/AleRax/ext/GeneRaxCore/ext/coraxlib'...
fatal: unable to access 'https://codeberg.org/Exelixis-Lab/coraxlib.git/': server certificate verification failed. CAfile: none CRLfile: none
```
You can skip the SSL certificate verification with:
```
git config --global http.sslVerify "false"
```
You can also have a look at this [post](https://forum.gitlab.com/t/gitlab-runner-server-certificate-verification-failed/59450/8) for a more detailed explanation of the problem.


To build the sources:
```
./install.sh
```

## Updating to the last stable version


```
./gitpull.sh # does a bit more than git pull
./install.sh
```

## Running

See the wiki (https://github.com/BenoitMorel/AleRax/wiki)

## Issues and questions

If you encounter any issue, please report it!!
For questions, issues or feedback, please post on the GeneRax google group: https://groups.google.com/g/generaxusers.
When reporting an issue, please send us at least the command line you ran, the logs file and the families file. The more information we get, the quicker we can solve the problems :-)

