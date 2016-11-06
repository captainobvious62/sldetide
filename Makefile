
# Build environment can be configured the following
# environment variables:
#   CC : Specify the C compiler to use
#   CFLAGS : Specify compiler options to use

CFLAGS += -I. -DPACKAGE_VERSION=\"1.0.1\"
LDFLAGS =
LDLIBS = -ltidal -ldali -lslink -lmseed -lm

all: sldetide msdetide

sldetide: sldetide.o
	$(CC) $(CFLAGS) -o $@ sldetide.o $(LDFLAGS) $(LDLIBS)

msdetide: msdetide.o
	$(CC) $(CFLAGS) -o $@ msdetide.o $(LDFLAGS) $(LDLIBS)

clean:
	rm -f sldetide.o sldetide msdetide.o msdetide

# Implicit rule for building object files
%.o: %.c
	$(CC) $(CFLAGS) -c $<

install:
	@echo
	@echo "No install target, copy the executable(s) yourself"
	@echo
