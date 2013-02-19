OBJDIR=objs
SRCDIR=src
DEPDIR=deps
CXX=mpic++

mpi:=$(shell which mpic++ 2>/dev/null)
ifeq ($(mpi),)
$(error "mpic++ not found - did you set your environment variables or load the module?")
endif

SRCS=\
     $(SRCDIR)/main.cpp \
     $(SRCDIR)/wsp.cpp \
     $(SRCDIR)/wsp_serial.cpp

OBJS=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS=$(patsubst $(SRCDIR)/%.cpp,$(DEPDIR)/%.d,$(SRCS))

CXXFLAGS+=-Wall -Wextra
LDFLAGS+=-lpthread -lmpi -lmpi_cxx

.PHONY: all lint clean pretty run

# all should come first in the file, so it is the default target!
all : wsp

run : wsp
	mpirun -np 4 wsp -i input/dist4

wsp: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

$(OBJS): | $(OBJDIR)
$(DEPS): | $(DEPDIR)
$(OBJDIR) $(DEPDIR):
	mkdir -p $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp cpplint.py Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -c -o $@
	python cpplint.py $< || true

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp Makefile
	@set -e; rm -f $@; \
        $(CXX) $(CXXFLAGS) -M $(CPPFLAGS) $< > $@.$$$$; \
        sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

-include $(DEPS)

clean:
	rm -rf $(OBJDIR) wsp

veryclean: clean
	rm -rf $(DEPDIR) cpplint.py scoreboard_token

lint : cpplint.py $(SRCS)
	python cpplint.py $(SRCS)

cpplint.py:
	wget http://google-styleguide.googlecode.com/svn/trunk/cpplint/cpplint.py
	chmod +x cpplint.py

pretty: $(SRCS)
	@# Check to make sure we have the proper version of indent.
	$(if $(findstring GNU,$(shell indent --version 2>/dev/null)),, $(error Install GNU indent))
	
	@echo "Running GNU indent on source files"
	@# Note - If we have any typedefs, we need to explicitly tell indent
	@# about our typedefs via
	@# 	-T type_t
	@indent \
	  -T solution_t \
	  --no-tabs \
	  --indent-level2 \
	  --continue-at-parentheses \
	  --line-length80 \
	  --honour-newlines  \
	  --blank-lines-after-procedures \
	  --break-before-boolean-operator \
	  --braces-on-if-line \
	  --braces-on-func-def-line \
	  --braces-on-struct-decl-line \
	  --cuddle-else \
	  --cuddle-do-while \
	  --case-brace-indentation0 \
	  --case-indentation0 \
	  --no-space-after-function-call-names \
	  --space-after-cast \
	  --no-space-after-parentheses \
	  --start-left-side-of-comments  \
	  --space-after-for \
	  --space-after-while \
	  --dont-break-procedure-type \
	  $(SRCS)
