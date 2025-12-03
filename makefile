# ==============================================================================
# VARIABLES DE CONFIGURATION
# ==============================================================================

# Compilateur
CXX = g++

# Répertoires
SRC_DIR = src
LIB_DIR = lib
BUILD_DIR = build
BIN_DIR = bin
INPUT_FILE = input/input.toml

# Nom de l'exécutable
EXEC = $(BIN_DIR)/ns_solver_2d

# --- Chemins d'Inclusion ---
# Adaptez ces chemins selon l'endroit où sont installées Eigen et la bibliothèque TOML
# Si Eigen est dans ${HOME}/libraries/eigen/, ce chemin est correct.
EIGEN_INCLUDE = ${HOME}/libraries/eigen/
TOML_INCLUDE = /net/netud/m/mdefaux/libraries/
# Options de compilation standard (C++17 pour la plupart des bibliothèques modernes de calcul)
CXXFLAGS = -std=c++17 -O2 -Wall -I$(LIB_DIR) -I$(EIGEN_INCLUDE) -I$(TOML_INCLUDE)

# Tous les fichiers source (y compris DataFile.cpp)
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Tous les fichiers objets correspondants (placés dans build/)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# ==============================================================================
# RÈGLES
# ==============================================================================

.PHONY: all debug clean run mkdirs

# Règle par défaut
all: release

# Compilation en mode Release (Optimisée)
release: CXXFLAGS += -DRELEASE
release: $(EXEC)
	@echo "--- Compilation RELEASE réussie. Exécutable: $(EXEC) ---"

# Compilation en mode Debug (inclut les symboles de débogage et pas d'optimisation)
debug: CXXFLAGS += -g -O0 -DDEBUG
debug: $(EXEC)
	@echo "--- Compilation DEBUG réussie. Exécutable: $(EXEC) ---"

# L'exécutable dépend de tous les objets
$(EXEC): $(OBJS) | mkdirs
	@echo "Liaison de l'exécutable..."
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation des objets
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | mkdirs
	@echo "Compilation de $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Création des dossiers si nécessaire
mkdirs: $(BUILD_DIR) $(BIN_DIR)
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# Nettoyage
clean:
	@echo "Nettoyage des fichiers résultats..."
	$(RM) -rf results/*.vtk
	$(RM) -rf results/*.dat

clean_all:
	@echo "Nettoyage des fichiers objets, de l'exécutable, et des dossiers..."
	$(RM) -rf $(BUILD_DIR) $(BIN_DIR)

# Compiler et exécuter avec le fichier d'entrée TOML
run: $(EXEC)
	@echo "Exécution de $(EXEC) avec $(INPUT_FILE)..."
	./$(EXEC) $(INPUT_FILE)