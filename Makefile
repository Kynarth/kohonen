.SILENT: clean

all: pdb dssp convert angle learning handling

handling:
	@echo "Handling kohonen data"
	@./scripts/traitement_matrice_kohonen.R
	@echo "Handling over"

learning:
	@echo "Learning on progress"
	@./scripts/matrice_de_kohonen.R
	@echo "Learning over"

angle: 
	@echo "Get angles"
	@python scripts/prepare_kohonen.py angle
	@echo "Complete recovery angles"

convert:
	@python scripts/prepare_kohonen.py convert -rv $(shell find dssp/ -type f | sort)

dssp: $(PDB_FILES)
	python scripts/prepare_kohonen.py dssp -v $(wildcard pdbs/list_pdb_ids/*.pdb)

pdb: list_pdb_ids.txt
	@python scripts/prepare_kohonen.py pdb -v $^

clean:
	@rm -rf pdbs dssp results/*.csv results/Kohonen/*.txt results/Kohonen/*.png\
		results/Kohonen/interm/csv/* results/Kohonen/interm/img/* \
		results/Kohonen/*.csv
