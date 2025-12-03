# Configuration de sortie (GIF Animé)
set terminal gif animate delay 10 size 800,800
set output 'simulation_vorticite.gif'

# Esthétique
set title "Simulation Navier-Stokes (Vorticité + Vitesse)"
set xlabel "X"
set ylabel "Y"
set size square
set xrange [0:1]
set yrange [0:1]

# Palette de couleurs (Bleu = Horaire, Blanc = Zéro, Rouge = Anti-horaire)
set palette defined (0 "blue", 1 "white", 2 "red")

# IMPORTANT : La vorticité peut être forte aux coins.
# Règle cette plage pour bien voir le contraste.
# Si c'est tout blanc, diminue la plage (ex: [-5:5]).
# Si c'est tout saturé, augmente la plage (ex: [-20:20]).
set cbrange [-10:10] 

# Configuration des vecteurs (Vitesse)
scale_factor = 0.05 

# Boucle sur les fichiers
fin = 4000 
step = 10 

do for [i=0:fin:step] {
    filename = sprintf('results/sol_%d.dat', i)
    
    # Vérification fichier
    check = system(sprintf("ls %s 2>/dev/null", filename))
    
    if (check ne "") {
        set title sprintf("Vorticité - Iteration : %d", i)
        
        # CHANGEMENT ICI : using 1:2:6 (colonne 6 est la vorticité)
        plot filename using 1:2:6 with image title "Vorticité", \
             filename using 1:2:($4*scale_factor):($5*scale_factor) with vectors head filled lc rgb "black" title "Vitesse"
    }
}