all: createIndex scala

.PHONY: createIndex scala libFMD libFMD-jar libsuffixtools

scala: libFMD-jar
	sbt stage

createIndex:
	$(MAKE) -C createIndex

libFMD:
	$(MAKE) -C libFMD

libFMD-jar:
	$(MAKE) -C libFMD jar-install

libsuffixtools:
	$(MAKE) -C libsuffixtools 
