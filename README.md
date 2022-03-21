# MeMC
A Monte-Carlo code for the simulation of fluctuating enclosed membranes. Such system can be
of relevance for the study of viruses, exosomes, e.t.c.

## Installation

If all the prequisites are satisfied, the installation is easy. The following
commands should do the job

'''bash
git clone https://github.com/vipinagrawal25/MeMC
cd MeMC
make
'''

## Prequisites

The code requires the following:

1)
Given a random configuration of particles distributed over a surface of sphere, the
code generates the minimum energy configuration. 

## Example

The surface config shown below::

![plot](./doc/figs/surf_mc_random.png)

After 60000 monte-carlo steps become::

![plot](./doc/figs/surf_mc_final.png)
