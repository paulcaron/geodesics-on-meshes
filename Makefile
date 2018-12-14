CACHE=.cache
all: verify

verify : build | $(CACHE)
	pmd -d src/ -cache $(CACHE) -f text -R ruleset.xml

$(CACHE):
	touch $@

build: clean
	javac -Xlint:all -cp lib/*:src/ src/* || exit 1

clean:
	rm -f src/*.class
