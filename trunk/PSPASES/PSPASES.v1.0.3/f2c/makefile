.SUFFIXES: .c .f

C_FILES = b_ax.c

F2C_OPTS = -I/usr/local/lam/include -c -w

all: $(C_FILES)

clean:
	rm $(C_FILES)

.f.c:
	f2c $(F2C_OPTS) $<
