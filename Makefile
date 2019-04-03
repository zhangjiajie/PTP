.PHONY: install develop clean

install:
	python3 setup.py install

develop: clean
	python3 setup.py develop --user

deb: clean
	touch tmp
	rm -rf tmp *.deb
	python3 setup.py --command-packages=stdeb.command bdist_deb

clean:
	python3 setup.py clean
	rm -rf build/ dist/ deb_dist/ *.egg-info/
	find . -name '*.pyc' -delete

dist:
	python3 setup.py bdist_egg

init:
	pip3 install -r requirements.txt
