

# AleRax  

AleRax is a parallel tool for species tree - gene tree inference and reconciliation under gene duplication, loss, and HGT. For each gene family, it takes as input a gene tree distribution (typically inferred with Bayesian inference tools such as MrBayes, PhyloBayes, etc.). AleRax can perform the following operations:
* Species tree inference
* Species tree rooting 
* Reconciled gene tree sampling 
* Model parameter estimation (e.g. DTL event probabilities) 
* Statistical test of different species tree hypotheses (you'll need to instal consel)
  
We are also working on the following features:
* Relative order of speciation event (relative dating) from HGT constraints
* Inference of highways of transfers (pairs of species that exchanged many genes via HGT)

When using AleRax, please cite: [https://academic.oup.com/bioinformatics/article/40/4/btae162/7633408](https://academic.oup.com/bioinformatics/article/40/4/btae162/7633408)

## Requirement

* A Linux or MacOS environnement
* gcc 5.0 or > 
* CMake 3.6 or >
* MPI (required if you want to use parallelization)

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

The generated executable is located here:
```
build/bin/alerax
```

To copy the executable to your PATH, such that you can call alerax from anywhere:
```
cd build
sudo make install
```

## Updating to the last stable version

In case we've made some changes since the last time you updated your repository, you can get and install the most recent version of AleRax with:

```
./gitpull.sh # does a bit more than git pull
./install.sh
```

Note that `git pull` does not update the submodules used by AleRax, hence the need for `./gitpull.sh`

## Running

See the wiki (https://github.com/BenoitMorel/AleRax/wiki)

## Issues and questions

If you encounter any issue, please report it!! I'm always happy to help.
For questions, issues or feedback, please post on the GeneRax (even if it's AleRax) google group: https://groups.google.com/g/generaxusers.
When reporting an issue, please send me at least the command line you ran, the logs file, and any file that might be relevant (e.g. the species tree file if AleRax failed to read the species tree...)


