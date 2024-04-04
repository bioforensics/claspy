## #===== development tasks =====#

## help:        print this help message and exit
help: Makefile
	@sed -n 's/^## //p' Makefile

## test:        run automated test suite
test:
	pytest --cov=claspy claspy

## format:      autoformat Python and Snakemake code
format:
	black --line-length=99 setup.py claspy/*.py claspy/tests/*.py

## style:       check code style
style:
	black --line-length=99 --check setup.py claspy/*.py claspy/tests/*.py

## hooks:       deploy git pre-commit hooks for development
hooks:
	echo "set -eo pipefail" > .git/hooks/pre-commit
	echo "make style" >> .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit
