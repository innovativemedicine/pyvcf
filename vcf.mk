ROOT := .
SCRIPTS := $(ROOT)/script
TESTS := $(ROOT)/test
3RDPARTY_SCRIPTS := $(ROOT)/3rdparty/script
PYTHON := python

export PYTHON ROOT

all: src/vcf/vcfparser.py

src/vcf/vcfparser.py: src/vcf/vcfparser.g
	$(3RDPARTY_SCRIPTS)/yapps2.py $<
	chmod +x $@
