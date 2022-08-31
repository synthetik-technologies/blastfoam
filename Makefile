# DESTDIR = /opt
prefix = /opt

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
	cp -r src $(DESTDIR)$(prefix)/src
	install --target-directory $(DESTDIR)$(prefix)/${BLAST_DIR}/platforms/${WM_OPTIONS}/etc -D \
		etc/bashrc
	install --target-directory $(DESTDIR)$(prefix)/${BLAST_DIR}/platforms/${WM_OPTIONS}/etc/codeTemplates -D \
		etc/codeTemplates/*

uninstall:
	rm -rf $(DESTDIR)$(prefix)
