all:
	cd code && $(MAKE) all
	cd tests && $(MAKE) all

clean:
	cd code && $(MAKE) clean
	cd tests && $(MAKE) clean



.PHONY: all clean
