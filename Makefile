all:
	make cython
	./test.py

clean:
	rm -rf __pycache__ *.c *.so build dist pyesg.egg-info

cython:
	python3 setup_cython.py build_ext --inplace
