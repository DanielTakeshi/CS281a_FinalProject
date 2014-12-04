# Daniel Seita

JFLAGS = -g
JC = javac

default: MainHarness.class

MainHarness.class: hmm/MainHarness.java
	$(JC) $(JFLAGS) hmm/MainHarness.java

clean:
	$(RM) *.class
	$(RM) hmm/*.class
	$(RM) utilities/*.class
