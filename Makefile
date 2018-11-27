all: build
	javac -cp lib/*:src/ src/*

build: clean

clean:
	rm -f src/*.class
