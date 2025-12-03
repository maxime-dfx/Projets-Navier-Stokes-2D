# ==========================================================
# Script Gnuplot : Animation de la Vorticité
# ==========================================================

# 1. Configuration Sortie GIF
set terminal gif animate delay 10 size 1200,600 optimize
set output 'Courbes/animation_vortex.gif'

# 2. Configuration Graphique 2D (Vue de dessus)
set view map
set size ratio -1  # Respecte le ratio X/Y (important pour ne pas écraser le cylindre)
set xrange [0:2.0]
set yrange [0:1.0]

# 3. Palette de couleurs (Bleu - Blanc - Rouge)
# Idéal pour la vorticité : Bleu (Rotation -), Blanc (0), Rouge (Rotation +)
set palette defined (-10 "blue", 0 "white", 10 "red")
set cbrange [-15:15] # Ajustez cette plage selon l'intensité de vos tourbillons !
set title "Vorticité (Curl U)"

# 4. Boucle sur les fichiers
# On suppose que vous avez des fichiers sol_0.dat, sol_20.dat ... sol_2000.dat
# Ajustez le '2000' et le pas '20' selon votre fréquence de sauvegarde

do for [i=0:10000:100] {
    file_name = sprintf('results/sol_%d.dat', i)
    set title sprintf("Temps : Allée de von Kármán (Iter %d)", i)
    
    # Colonne 1:x, 2:y, 6:omega (Vorticité)
    # 'with image' crée la heatmap pixelisée
    plot file_name using 1:2:6 with image notitle
}

print "Animation générée : Courbes/animation_vortex.gif"