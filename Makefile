all: createIndex speedMap

.PHONY: libsuffixtools libFMD libFMD-jar createIndex speedMap scala

scala: libFMD-jar
	sbt stage

speedMap:
	$(MAKE) -C speedMap

createIndex:
	$(MAKE) -C createIndex

libFMD:
	$(MAKE) -C libFMD

libFMD-jar:
	$(MAKE) -C libFMD jar-install

libsuffixtools:
	$(MAKE) -C libsuffixtools 
