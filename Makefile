all:
	python3 setup.py install

clean:
	rm -rf __pycache__ *.c *.so build dist pyesg.egg-info
