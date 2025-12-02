CC=gcc
CFLAGS=-I.

SRC_FILES = bin_to_txt.c txt_to_bin.c read_config.c
OBJ_FILES = $(SRC_FILES:.c=.o)
BIN_FILES = bin_to_txt txt_to_bin 
all: $(BIN_FILES)


%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

bin_to_txt: bin_to_txt.o read_config.o
	$(CC) -o $@ $^ $(CFLAGS)

txt_to_bin: txt_to_bin.o 
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f $(OBJ_FILES) $(BIN_FILES)
