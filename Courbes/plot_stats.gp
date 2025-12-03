# ==========================================================
# Script Gnuplot : Statistiques Allée de von Kármán
# ==========================================================

# 1. Sortie Image
set terminal pngcairo size 1000,600 enhanced font 'Arial,12'
set output 'Courbes/simulation_stats.png'

# 2. Titres et Légendes
set title "Simulation : Allée de von Kármán (Re ~ 120)"
set xlabel "Temps (s)"
set grid

# 3. Configuration des axes
set ylabel "Énergie Cinétique (J)" textcolor rgb "blue"
set ytics textcolor rgb "blue" nomirror

set y2label "Vitesse Max (m/s)" textcolor rgb "red"
set y2tics textcolor rgb "red"
set y2range [0:3.0] # On fixe un peu l'échelle pour voir les pics

# 4. Fichier de données
FILE = 'results/history.dat'

# 5. Tracer
# On s'attend à voir des oscillations régulières après t=3s
plot FILE using 1:2 axes x1y1 with lines lw 2 lc rgb "blue" title "Énergie Cinétique", \
     FILE using 1:3 axes x1y2 with lines lw 1 lc rgb "red" dashtype 2 title "Vitesse Max"

print "Statistiques générées : results/simulation_stats.png"