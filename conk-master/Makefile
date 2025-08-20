all : conk/conk.pyx setup.py README.md
	python3 setup.py sdist bdist_wheel
	pip3 install --user dist/conk*whl

clean :
	rm -rf build/ dist/ conk.egg-info/ conk/conk.c conk/__pycache__
