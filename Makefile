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
	install -D etc/bashrc /etc$(DESTDIR)$(prefix)/bashrc

uninstall:
	rm -rf $(DESTDIR)$(prefix)
	rm -rf /etc$(DESTDIR)$(prefix)