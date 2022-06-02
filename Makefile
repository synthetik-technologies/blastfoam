DESTDIR = /opt
prefix = /blastfoam

build: SHELL:=bash
build:
	+ source /opt/openfoam9/etc/bashrc && \
	source etc/bashrc && \
	./Allwmake -j -s

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
	# * find better install location later
	install --target-directory $(DESTDIR)$(prefix)/etc -D \
		etc/bashrc
	install --target-directory $(DESTDIR)$(prefix)/etc/codeTemplates -D \
		etc/codeTemplates/*

uninstall:
	rm -rf $(DESTDIR)$(prefix)
