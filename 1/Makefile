CC=gcc
.PHONY: calculate jacobi gs sor clean purge

# Compila i enllaça els tres programes, els executa, i emmagatzema els
# resultats en els fitxers jacobi.txt, gs.txt i sor.txt.
calculate: jacobi gs sor
	./jacobi.exe > jacobi.txt
	./gs.exe > gs.txt
	./sor.exe > sor.txt
	make clean

# Compila i enllaça jacobi.c.
jacobi: jacobi.o common.o
	$(CC) -o jacobi.exe jacobi.c

# Compila i enllaça gs.c.
gs: gs.o common.o
	$(CC) -o gs.exe gs.c

# Compila i enllaça sor.c.
sor: sor.o common.o
	$(CC) -o sor.exe sor.c

# Esborra tots els fitxers generats excepte els *.txt.
clean:
	rm -f *.exe *.o

# Esborra tots els fitxers generats incloent els *.txt.
purge:
	rm -f *.exe *.o jacobi.txt gs.txt sor.txt

%.o: %.c
	$(CC) -c -o $@ $<
