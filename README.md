This covers basically everything that's important to know about this code.

(c) November/December 2014 by Daniel Seita


***The Plan***

Here's my goal. I want to develop a technique to detect anomalous events that are generated via a
HMM. In other words, I'd like to implement what's described in the reference paper (1). I believe
what this means is that I will need to design my own HMMs and run them to generate some synthetic
data. This data can be either benign (normal) or dangerous (anomalous) and the goal is to use my
technique to detect with high confidence when the anomaly occurs. I believe I will need to mix up
the generated data somehow. As an example, consider the nuclear weapons experiment in the reference
paper (1). 

FIRST STEP ... is to figure out how to code up the HMMs so that I can generate a sequence of
observations. From the paper, it seems like I will need an HMM for each adversarial activity I have
(one is enough to start out) and ... I will need something else for the "noise". What that is, is
not clear, so I will just have to decide for myself.

THEN ... once I have that set up, I can implement Page's algorithm. But the point is to get the
sequence of observations set up now.

***Code Outline***

I think I should put everything in one package. Then there will be several programs in there,
described below. Most integration will be TBD ... basically I'll fix them up as they come along.
Sloppy, I know. Here are some ideas right now:

(1) class for nodes/states (the hidden ones)
(2) class for transactions (which use nodes/states)
(3) class for observations
(4) class for Page's test?
(5) class to actually run the HMM
(6) class for features
(7) class for various utilities


***References***

(1) Anomaly Detection via Feature-Aided Tracking and Hidden Markov Models. This is the main paper
that I want to be basing my project on.

(x) HMM-Based Reconstruction of Unreliable Spectrographic Data for Noise Robust Speech Recognition.
This is a nice application paper but I don't think I'll be using much of their stuff. Perhaps I
should cite this in a "Related Work" section.

(x) Robust Sparse Regression under Adversarial Corruption. This was the one Ben recommended to me
but I think I'm moving away from that idea. And besides, I'm using HMMs, which don't really need the
"trimmed inner product" that they use, I think.

(x) (Forgot name) There was a survey paper from the University of Minnesota that talked about how
HMMs are used in "anomaly detection." That might be useful to learn about background work.
