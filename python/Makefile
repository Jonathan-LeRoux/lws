default: build

############### FILES #################
FILES=lwslib lws.bycython.cpp setup.py MANIFEST.in README.rst

lwslib: ../lwslib
	cp -R ../lwslib lwslib

lws.bycython.cpp: lws.pyx lws_functions.pxd
	cython --cplus lws.pyx -o lws.bycython.cpp
#######################################

############## COMMANDS ###############
build: ${FILES}
	python setup.py build_ext -i
	touch build

sdist: ${FILES}
	python setup.py sdist
	touch dist

publish: clean sdist
	twine upload dist/*

clean:
	rm -rf lwslib dist lws.egg-info build
	rm -f lws.cpp *.bycython.* lws.*.so
#######################################
