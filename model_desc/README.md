
# Model description

## Fully-connected SIR

The population is composed of <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> indistinguishable individuals,
each interacting equally with all the others.  An individual can be in
one of three states S (susceptible), I (infected) or R (recovered). 
The dynamics is defined by the transition rates (per individual)

<p align="center"><img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/571f946a21eb2b3a2de077f5f1abc84d.svg?invert_in_darkmode" align=middle width=130.0391466pt height=63.5745132pt/></p>

The rate constants <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.16555099999999pt height=22.831056599999986pt/> and <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode" align=middle width=9.423880949999988pt height=14.15524440000002pt/> are related to the infection
time <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/7eb44e3c5f2ccb05ae52bafea039a0b4.svg?invert_in_darkmode" align=middle width=21.91224089999999pt height=20.221802699999984pt/> and <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode" align=middle width=19.034022149999988pt height=22.465723500000017pt/> through

$<img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/1773c59dba0725c6f84607dfa23930a5.svg?invert_in_darkmode" align=middle width=160.2571509pt height=27.77565449999998pt/>

The [./sir_par.dat][parameter file] gives <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode" align=middle width=19.034022149999988pt height=22.465723500000017pt/> and <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/7eb44e3c5f2ccb05ae52bafea039a0b4.svg?invert_in_darkmode" align=middle width=21.91224089999999pt height=20.221802699999984pt/> as
well as the initial conditions <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/7ea1d381135a7e0a3cd0eefafea7c973.svg?invert_in_darkmode" align=middle width=16.632471899999988pt height=22.465723500000017pt/> and <img src="https://rawgit.com/in	git@github.com:tgrigera/COVIDm/None/svgs/88fbd05154e7d6a65883f20e1b18a817.svg?invert_in_darkmode" align=middle width=13.77859724999999pt height=22.465723500000017pt/>, namely the initial
fraction of susceptible and infected.

