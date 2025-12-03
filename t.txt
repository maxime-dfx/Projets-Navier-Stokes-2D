# ==========================================================
# Script Gnuplot pour Navier-Stokes 2D
# ==========================================================

# 1. Configuration de la sortie (Image PNG propre)
set terminal pngcairo size 1000,600 enhanced font 'Arial,12'
set output 'results/simulation_stats.png'

# Optionnel : Affichage interactif (décommentez pour voir à l'écran)
# set terminal wxt size 1000,600 persist

# 2. Titres et Grille
set title "Simulation Navier-Stokes : Décroissance de Taylor-Green"
set xlabel "Temps (s)"
set grid

# 3. Configuration du double axe Y
# Axe de Gauche (Bleu) : Énergie Cinétique
set ylabel "Énergie Cinétique Totale" textcolor rgb "blue"
set ytics textcolor rgb "blue" nomirror

# Axe de Droite (Rouge) : Vitesse Max
set y2label "Vitesse Maximum (Magnitude)" textcolor rgb "red"
set y2tics textcolor rgb "red"

# 4. Échelle Logarithmique (Optionnel)
# Taylor-Green est censé décroître exponentiellement.
# Décommentez la ligne suivante pour vérifier si vous obtenez une droite :
# set logscale y

# 5. Chemin du fichier de données
# ATTENTION : Vérifiez que c'est bien le dossier défini dans votre .toml
FILE = 'results/history.dat'

# 6. Tracé
# using 1:2 = t vs Ek (Axe x1y1)
# using 1:3 = t vs Vmax (Axe x1y2)
plot FILE using 1:2 axes x1y1 with lines lw 3 lc rgb "blue" title "Énergie Cinétique", \
     FILE using 1:3 axes x1y2 with lines lw 2 lc rgb "red" dashtype 2 title "Vitesse Max"

print "Graphique généré : results/simulation_stats.png"