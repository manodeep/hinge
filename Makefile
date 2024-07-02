include common.mk

DIRS = haloparentfinder orphanfixer mergertree
BUILDDIRS = $(DIRS:%=build-%)
CLEANDIRS = $(DIRS:%=clean-%)

all: $(BUILDDIRS)

tar:
	git archive --format=tar.gz -o $(DISTNAME).$(MAJOR).0.$(MINOR).no_version_control.tar.gz --prefix=hinge/ master

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
