include common.mk

DIRS = haloparentfinder orphanfixer mergertree
BUILDDIRS = $(DIRS:%=build-%)
CLEANDIRS = $(DIRS:%=clean-%)

all: $(BUILDDIRS)

tar:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).no_version_control.tar.gz -X ".hg*"

dist:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz 

.PHONY: clean celna clena celan $(DIRS) $(BUILDDIRS) $(CLEANDIRS) 

$(DIRS): $(BUILDDIRS)
$(BUILDDIRS):
	$(MAKE) -C $(@:build-%=%)

clean: $(CLEANDIRS)
	$(RM) $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz

$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

clena: clean
celan: clean
celna: clean


