CC     = mpicc
CFLAGS = -std=c99 -Wall -O3
LIBS   = -lmpi -lm

wave: wave.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	rm wave

