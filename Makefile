CFLAGS = -Wall
SACLIB = /opt/SAC/lib

BIN = ./bin

all: sac2col sacch saclh sacmax sacstack clean

sac2col: sac2col.o sacio.o
	$(CC) -o $(BIN)/$@ $^

sacch: sacch.o sacio.o datetime.o
	$(CC) -o $(BIN)/$@ $^ -lm

saclh: saclh.o sacio.o
	$(CC) -o $(BIN)/$@ $^

sacmax: sacmax.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -lm

sacstack: sacstack.o sacio.o
	$(CC) -o $(BIN)/$@ $^ -lm -L$(SACLIB) -lsac -lsacio

clean:
	rm *.o
