---------------------------------------------------------------------------------------------------

HMM-Detector is a program that can detect changes in HMM distributions.
Copyright (C) 2014 Daniel Seita

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not,
see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------------------------------

*Some short-term goals*

(1) Continue to do research to see if anyone's done something like this.
(2) Study the sequential likelihood ratio test and the cumulative sum test.
(3) Decide on what experiments to run.

***The Plan***

My goal in doing this project is to analyze statistical techniques that can detect a sudden change
in a distribution based on a set of Hidden Markov Models (HMMs).  Here is the ideal setting: we have
k different HMMs, each with their own distribution that are pairwise "roughly similar" yet
"sufficiently different" with each other. (I know this sounds vague, but part of my goal will be to
research how much I can make them "similar" or "different.) We will "run" these HMMs in some order,
probably one complete run, then another complete run. As an example, imagine that we have two HMMs
that we'll call A and B. We will output a series of observations based on going through A, from its
start state to its end state. Then, right after that last observation, we immediately switch over to
B, so the full cycle is: start at A's start state, proceed through the HMM, go to A's end state,
then go to B's start state, proceed through the new HMM, and then go to B's end state, at which time
we terminate the run. The most important goal will be to determine with high confidence when the
distribution changes from being based on A to being based on B. To be specific, let's index our
observations by

a1 a2 ... am b(1+m) b(2+m) ... b(n+m)

Where a1 = the start state for HMM A, etc. Then our statistical method might return an index k and
the goal will be to have k as close to 1+m as possible. We might enforce more constraints if k <
1+m, i.e., we could do this by minimizing Bayes' risk and assign costs appropriately.

Here are some interesting research questions based off of this setting:

1. (More about what I should be doing) What should the statistical method do? What should it try to
learn or minimize (e.g., Bayes' risk)? Should it return a probability of whether we're at the new
HMM after each observation? Of course, determining appropriate costs is a related issue. I'll
probably need to fix one way of doing this (at most two) and assume it for all experiments.

2. What kind of statistical test should we be running to detect the change in distribution based on
this setting? The reference paper uses Page's test. What about if we use a different technique, say
based on data-mining, maximum likelihood? They are not completely clear on how they implemented that
alternative.

3. When using Page's test, they took the conditional likelihood ratio between the reference/null
distribution and the target/"anomaly" distribution. Why not just use the joint distribution? Is this
kind of question similar to what Andrew Ng did in his Generative vs. Discriminative paper?

4. What about if we modify the priors of the HMM? How much strength does the prior distribution have
on when we should declare a change in distributions?

5. How much can we effectively vary the h parameter? (In the reference paper, if the cumulative sum
statistic was greater than h, a "change in distribution" was declared and the method terminated.)
Would it be worth to develop a tuning method?

6. In Appendix III, they discuss several speed-ups related to skipping over null observations. But I
don't understand how to apply that here, because we'd have to *know* that the observations are null
(i.e., come from the non-anomalous HMM), but we don't know for sure because everything is run
through *hidden* states. So I am curious as to whether there are other speedups that are more
effective to implement,

7. (*MAJOR*) How should I run experiments? What kind of distributions do I create? How many should I
create, and in what order do I run them? What are the best factors to vary or keep the same? How can
I know if the test is working, because two different distributions vary wildly in how much they
differ, so a weaker test applied to two very different distributions can detect a change quickly. I
will need to think carefully about this.

My original goal in this project was to develop a technique to detect anomalous events that are
generated via a HMM, with the help of "feature tracking.  In my opinion, modeling feature tracking
is not that important here because that is contained in the transition/emission probabilities.  If I
need to "pretend" that a new HMM is similar to a reference one but different in a few aspects (i.e.,
"features"), I can set that HMM to have the same number of states and just tweak the
transition/emission probabilities from the reference's stats. And also, I don't see a need to model
an HMM by using observations based on two states, but it's something worth exploring. I had so many
questions based on reading the reference paper that I think it's fair I get to test out some of the
choices they made. I will also need to fill in some details not explained in the paper.


***Code Outline***

Right now this will be in perpetual draft state until I get something decent going. But in the
meantime, here's what I think:

(1) I'll need a program that can compute forward probabilities (or other important computations,
such as the scaled forward they mention in the paper) for generic HMMs. This should probably be the
program that does Page's test.

(2) I'll need a main harness that can do the whole pipeline. In particular, it will be running the
experiments (or I can put the experiments in a separate class if it gets too "bulky").

(3) I'll probably want a generic HMM interface so that I can make one HMM class (that implements
this interface) for each HMM I want to build.


***References***

Goal: have many citations from 2009-2014 as possible.

(1) Anomaly Detection via Feature-Aided Tracking and Hidden Markov Models. This is the main paper
that I want to be basing my project on. I should also look at David Willet's more recent work that
he's done, though his website isn't up to date. Ideally, I should say that much of this work is
unexplored.

(2) That popular 1989 reference paper on HMMs: "A tutorial on hidden Markov models and selected
applications in speech recognition" by L. R. Rabiner.

(3) HMM-Based Reconstruction of Unreliable Spectrographic Data for Noise Robust Speech Recognition.
This is a nice application paper but I don't think I'll be using much of their stuff. Perhaps I
should cite this in a "Related Work" section.

(4) Robust Sparse Regression under Adversarial Corruption. This was the one Ben recommended to me
but I think I'm moving away from that idea. And besides, I'm using HMMs, which don't really need the
"trimmed inner product" that they use, I think.

(5) (Forgot name) There was a survey paper from the University of Minnesota that talked about how
HMMs are used in "anomaly detection." That might be useful to learn about background work.
Referencing a couple of survey papers probably won't hurt.
