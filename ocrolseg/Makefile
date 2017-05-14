CXXFLAGS=-fPIC -g -O4

all: libocrolseg.a ocrolseg.py

install: all
	sudo python setup.py install

ocrolseg.py: ocrolseg.i
	python setup.py build_ext

libocrolseg.a: ocrolseg.o
	ar cr $@ $^

clean:
	rm -rf build *.o *.so *.a *_wrap.cpp ocrolseg.py
