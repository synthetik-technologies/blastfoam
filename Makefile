DESTDIR = /opt
prefix = /blastfoam

build: SHELL:=bash
build:
	+ source /opt/openfoam9/etc/bashrc && \
	source etc/bashrc && \
	./Allwmake -j

clean: SHELL:=bash
clean:
	+ source /opt/openfoam9/etc/bashrc && \
	source etc/bashrc  && \
	./Allwclean

install:
	install --target-directory $(DESTDIR)$(prefix)/bin -D \
		bin/*
	install --target-directory $(DESTDIR)$(prefix)/lib -D \
		lib/*
	install --target-directory $(DESTDIR)$(prefix)/etc -D \
		etc/* #* temporary, find better placement later

uninstall:
	rm -rf $(DESTDIR)$(prefix)
