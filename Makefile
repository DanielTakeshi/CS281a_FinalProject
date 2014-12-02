# Daniel Seita

JFLAGS = -g
JC = javac

default: RunHMM.class

RunHMM.class: hmm/RunHMM.java
	$(JC) $(JFLAGS) hmm/RunHMM.java

clean:
	$(RM) *.class
	$(RM) hmm/*.class
	$(RM) utilities/*.class
