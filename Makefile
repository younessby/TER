# Nom du compilateur Fortran
FC = gfortran

# Flags de compilation, incluant -g pour les informations de débogage
FFLAGS = -O2 -g

# Nom du programme exécutable
PROGRAM = programme

# Liste des modules
MODULES = points.o arretes.o triangles.o triangulation.o

# Liste des objets
OBJECTS = $(MODULES) programme.o

# Règle par défaut
all: $(PROGRAM)

# Règle pour créer le programme exécutable
$(PROGRAM): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

# Règle pour les modules et le programme principal
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Nettoyage
clean:
	rm -f $(PROGRAM) $(OBJECTS)

# Dépendances
points.o: points.f90
arretes.o: arretes.f90 points.o
triangles.o: triangles.f90 arretes.o points.o
triangulation.o: triangulation.f90 points.o arretes.o triangles.o
programme.o: programme.f90 triangulation.o points.o arretes.o triangles.o
