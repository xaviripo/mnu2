CC=gcc
.PHONY: main clean purge

# Compila i enllaça el programa i l'executa.
calculate: main
	./main.exe
	make clean

# Compila i enllaça main.c.
main: main.o
	$(CC) -o main.exe main.c -lm

# Esborra tots els fitxers generats excepte el *.txt.
clean:
	rm -f *.exe *.o

# Esborra tots els fitxers generats incloent el *.txt.
purge:
	rm -f *.exe *.o output.txt

# Compila els fitxers font a fitxers objecte
%.o: %.c
	$(CC) -c -o $@ $<
