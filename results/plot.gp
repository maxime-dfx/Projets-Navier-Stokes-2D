# Configuration de sortie (GIF Animé)
set terminal gif animate delay 10 size 800,800
set output 'simulation.gif'

# Esthétique
set title "Simulation Navier-Stokes (Pression + Vitesse)"
set xlabel "X"
set ylabel "Y"
set size square
set xrange [0:1]
set yrange [0:1]

# Palette de couleurs (Bleu -> Blanc -> Rouge) pour la Pression
set palette defined (0 "blue", 1 "white", 2 "red")
set cbrange [0:1] # ADAPTEZ CES VALEURS selon vos résultats (grep) !

# Configuration des vecteurs (Vitesse)
# scale : taille des flèches (à réduire si trop grosses)
scale_factor = 0.05 

# Boucle sur les fichiers .dat
# Adaptez 'fin' selon le nombre de fichiers générés
fin = 4000 # Nombre approximatif d'itérations sauvegardées
step = 10 # Votre 'save_freq' dans le main.cpp

do for [i=0:fin:step] {
    filename = sprintf('results/sol_%d.dat', i)
    
    # Vérifie si le fichier existe (commande Linux)
    check = system(sprintf("ls %s 2>/dev/null", filename))
    
    if (check ne "") {
        set title sprintf("Temps iteration : %d", i)
        
        # Le plot magique :
        # 1. 'with image' utilise la colonne 3 (Pression) pour la couleur
        # 2. 'with vectors' utilise colonnes 1,2 (XY) et 4,5 (UV) pour les flèches
        
        plot filename using 1:2:3 with image title "Pression", \
             filename using 1:2:($4*scale_factor):($5*scale_factor) with vectors head filled lc rgb "black" title "Vitesse"
    }
}