# You can include this file in your project's makefile like so:
#
# VPATH = path/to/pyvcf
# include path/to/pyvcf/vcf.mk
# something.py: src/vcf/vcfparser.py

PYTHON := python
YAPPS := 3rdparty/script/yapps2.py

export PYTHON ROOT

all: src/vcf/vcfparser.py

src/vcf/vcfparser.py: src/vcf/vcfparser.g $(YAPPS)
	$(word 2,$^) $<
	$(eval VCFPARSER := $$(patsubst %.g,%.py,$$<))
	chmod +x $(VCFPARSER)
