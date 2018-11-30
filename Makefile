#Makefile
#Reference  https://qiita.com/szkny/items/07b4c93702b3f8090fa3 


SUFFIX   = .cpp
COMPILER = g++
CFLAGS   = -Ofast

SRCDIR   = ./source
INCLUDE  = ./include 
EXEDIR   = ./bin

SOURCE   = $(wildcard $(SRCDIR)/*$(SUFFIX))
OBJECTS  = $(notdir $(SOURCE:$(SUFFIX)=.o))
TARGETS  = $(notdir $(basename $(SOURCE)))

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

define MAKEALL
$(1): $(1).o
	$(COMPILER) -I $(INCLUDE) $(CFLAGS)  -o $@ $(EXEDIR)/$(1) $(1).o $(ROOTLIBS)  
	@$(RM) $(1).o
$(1).o:
	$(COMPILER) $(ROOTFLAGS)  -I $(INCLUDE) $(CFLAGS) -c $(SRCDIR)/$(1)$(SUFFIX)
endef

.PHONY: all
all: $(TARGETS)
$(foreach VAR,$(TARGETS),$(eval $(call MAKEALL,$(VAR))))

#make clean
.PHONY: clean
clean:
	$(RM) $(EXEDIR)/* *.o
