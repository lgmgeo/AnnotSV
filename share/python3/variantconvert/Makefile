VERSION=$$(grep '__version__ =' src/variantconvert/__init__.py | cut -d '"' -f2)

# use run like this:
# $(make run) -args -for --variantconvert
# otherwise makefile thinks args are for itself
run:
	@echo python src/variantconvert/__main__.py

#other routines don't need args and can be used normally
build:
	python -m build

pypi:
	python -m twine upload --repository testpypi dist/variantconvert-$(VERSION)*;

install:
	python -m build
	pip install dist/variantconvert-$(VERSION).tar.gz
	variantconvert init

#a current issue I'm trying to fix
debug:
	python src/variantconvert/__main__.py convert -i tests/data/DECON.results_all.AnnotSV.tsv -o decon_annotsv_test.vcf -c src/variantconvert/configs/hg19/annotsv3_from_vcf.json -v debug